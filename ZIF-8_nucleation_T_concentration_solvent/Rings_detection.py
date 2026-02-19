import numpy as np
import networkx as nx

# ====================================
zn_type = 1          # Atom type for Zn
n1_type = 3          # Atom type for N (assuming ligands contain N)
max_zn_n_distance = 3.0  # Maximum Zn-N distance to consider bonding 

traj_file = "Trajectory.lammpstrj"
count_file = "count_rings_based-time.dat"
ring_file = "ring_atoms_indexes_based-time.dat"
matrix_file = "ring_matrix_based-time.dat"

# Time parameters
deltaT = 25        # ps
offsetT = 0.0      # ps
# ====================================

def read_frame(file):
    global box_bounds
    atoms = []
    while True:
        line = file.readline()
        if not line:
            return None
        if line.startswith("ITEM: TIMESTEP"):
            timestep = int(file.readline().strip())
        elif line.startswith("ITEM: NUMBER OF ATOMS"):
            n_atoms = int(file.readline().strip())
        elif line.startswith("ITEM: BOX BOUNDS"):
            box_bounds = []
            for _ in range(3):
                bounds = list(map(float, file.readline().split()))
                box_bounds.append((bounds[0], bounds[1]))
        elif line.startswith("ITEM: ATOMS"):
            for _ in range(n_atoms):
                parts = file.readline().split()
                atom_id = int(parts[0])
                atom_type = int(parts[1])
                x, y, z = map(float, parts[3:6])
                atoms.append((atom_id, atom_type, np.array([x, y, z])))
            return timestep, atoms

def minimum_image_distance(pos1, pos2, box_bounds):
    delta = pos1 - pos2
    box_lengths = np.array([high - low for (low, high) in box_bounds])
    for i in range(3):
        if delta[i] > 0.5 * box_lengths[i]:
            delta[i] -= box_lengths[i]
        elif delta[i] < -0.5 * box_lengths[i]:
            delta[i] += box_lengths[i]
    return np.linalg.norm(delta)

def extract_zn_and_ligands(atoms, box_bounds):
    zn_atoms = []
    n_atoms = []
    
    for atom_id, atom_type, pos in atoms:
        if atom_type == zn_type:
            zn_atoms.append((atom_id, pos))
        elif atom_type == n1_type:  # 
            n_atoms.append((atom_id, pos))
    
    # Find N pairs that belong to the same ligand 
    ligands = []
    used_n = set()
    
    for i, (n1_id, n1_pos) in enumerate(n_atoms):
        if n1_id in used_n:
            continue
        # Find the closest N to n1 
        min_dist = float('inf')
        closest_n2 = None
        
        for j, (n2_id, n2_pos) in enumerate(n_atoms[i+1:], i+1):
            if n2_id in used_n:
                continue
            dist = minimum_image_distance(n1_pos, n2_pos, box_bounds)
            if dist < min_dist and dist < 2.5:  # 
                min_dist = dist
                closest_n2 = (n2_id, n2_pos)
        
        if closest_n2:
            ligands.append({
                'n1_id': n1_id,
                'n1_pos': n1_pos,
                'n2_id': closest_n2[0],
                'n2_pos': closest_n2[1]
            })
            used_n.add(n1_id)
            used_n.add(closest_n2[0])
    
    return zn_atoms, ligands

def build_zn_zn_graph(zn_atoms, ligands, box_bounds):
    G = nx.Graph()
    G.add_nodes_from(zn_id for zn_id, _ in zn_atoms)

    for ligand in ligands:
        n1_pos = ligand['n1_pos']
        n2_pos = ligand['n2_pos']

        # Find closest Zn to n1
        min_dist1 = float('inf')
        closest_zn1 = None
        for zn_id, zn_pos in zn_atoms:
            dist = minimum_image_distance(n1_pos, zn_pos, box_bounds)
            if dist < min_dist1 and dist < max_zn_n_distance:
                min_dist1 = dist
                closest_zn1 = zn_id

        # Find closest Zn to n2
        min_dist2 = float('inf')
        closest_zn2 = None
        for zn_id, zn_pos in zn_atoms:
            dist = minimum_image_distance(n2_pos, zn_pos, box_bounds)
            if dist < min_dist2 and dist < max_zn_n_distance:
                min_dist2 = dist
                closest_zn2 = zn_id

        if closest_zn1 is not None and closest_zn2 is not None and closest_zn1 != closest_zn2:
            G.add_edge(closest_zn1, closest_zn2)

    return G


# Initialize output files
with open(count_file, 'w') as f:
    f.write("# Time(ps) 3R 4R 5R 6R 7R 8R\n")
with open(ring_file, 'w') as f:
    f.write("# Ring sizes and atom indices\n")
with open(matrix_file, 'w') as f:
    f.write("# Time(ps) 3R 4R 5R 6R 7R 8R\n")

# === Main loop ===
with open(traj_file) as f, open(count_file, "w") as count_out, open(ring_file, "w") as ring_out:
    frame = 0
    while True:
        result = read_frame(f)
        if result is None:
            break
            
        timestep, atoms = result
        current_time = frame * deltaT + offsetT


        zn_atoms, ligands = extract_zn_and_ligands(atoms, box_bounds)
        G = build_zn_zn_graph(zn_atoms, ligands, box_bounds)

        all_rings = nx.minimum_cycle_basis(G)

        valid_rings = [ring for ring in all_rings if 3 <= len(ring) <= 8]

        counts = {3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0}
        for ring in valid_rings:
            counts[len(ring)] += 1


        # Using Fortran-style write format for time
        count_out.write(f"{current_time:.6f} {counts[3]} {counts[4]} {counts[5]} {counts[6]} {counts[7]} {counts[8]}\n")

        ring_out.write(f"Time = {current_time:.6f} ps\n")
        for ring in valid_rings:
            ring_out.write(f"  {len(ring)}-member ring: " + " ".join(map(str, ring)) + "\n")

        # Write to matrix file using Fortran-style format
        with open(matrix_file, "a") as mfile:
            mfile.write(f"{current_time:.6f} " + " ".join(str(counts[r]) for r in range(3, 9)) + "\n")

        frame += 1

print(f"Done! Ring counts in '{count_file}', ring atoms in '{ring_file}', matrix in '{matrix_file}'")
