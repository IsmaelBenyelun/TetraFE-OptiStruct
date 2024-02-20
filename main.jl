"""Script to mesh a domain with tetrahedra and run static case via OptiStruct
"""

push!(LOAD_PATH, pwd())
using Mesh
using ParaView
using MyStructs
using WriteUtils
using OptiStruct


## Define global properties
DOF_PER_NODE = 3
E = 200.0e3 # MPa
ν = 0.33
F = 1.e3 # N
MATERIAL = Isotropic(E, ν)


### MESH -- cubic domain with tetrahedral elements
ELEM_TYPE = Tetra4()
elems_per_dim = [5, 5, 5]
sizes = [100., 100., 100.] # X, Y, Z

n_nodes, n_elems, v_nodes, v_mesh, mechanical_mesh = MeshBuildTetra(sizes, elems_per_dim, E, ν)
n_dofs = DOF_PER_NODE * n_nodes

# Save in .vtu format
if !isdir("paraview-output")
    mkdir("paraview-output")
end
PrintParaview(v_nodes, v_mesh, "./paraview-output/mesh.vtu")


## BOUNDARY CONDITIONS -- Hard-coded on purpose
force_nodes = [2, 4, 6, 8]
# force_nodes = [89, 95, 125, 131] .+ 1
force_dofs = force_nodes .* 3

restricted_nodes = [1, 3, 5, 7]
# restricted_nodes = [0, 6, 12, 18, 24, 30, 36, 66, 72, 102, 108, 138, 144, 174, 180, 186, 192, 198, 204, 210] .+ 1
restricted_dofs = reduce(vcat, [[i*3 - 2, i*3 - 1, 3*i] for i in restricted_nodes])

u = zeros(n_dofs)
f = zeros(n_dofs)
f[force_dofs] .= F


## OPTISTRUCT INTERFACE
SOLVER = "OPTISTRUCT"
solver_path = "D:\\Altair\\2021.2\\hwsolvers\\scripts\\optistruct.bat" # Absolute path to OptiStruct's executable file (optistruct.bat)
# Create .fem file
fileName = write_input_file(SOLVER, ELEM_TYPE, MATERIAL, v_nodes, v_mesh, mechanical_mesh, DOF_PER_NODE, u, f, restricted_dofs, force_dofs)
# Run .exe
u, ε_el, σ_el, ip_ID = optistruct(fileName, solver_path, DOF_PER_NODE)


## SAVE PARAVIEW FILE
PrintParaview(v_nodes, v_mesh, ε_el, σ_el, u, "./paraview-output/3d_tetra.vtu");
