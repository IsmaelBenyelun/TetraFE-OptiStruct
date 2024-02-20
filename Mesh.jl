module Mesh

using DelimitedFiles
export MeshBuildTetra


function MeshBuildTetra(v_size::Array{Float64, 1}, v_divs::Array{Int64,1}, E::Float64, ν::Float64)
    """ SIZE AND NUMBER OF DIVISIONS """
    # Unpack
    sizex, sizey, sizez = v_size
    ndivx, ndivy, ndivz = v_divs

    # Compute gap
    gapx = sizex / ndivx
    gapy = sizey / ndivy
    gapz = sizez / ndivz

    # Nodes in each directions
    nodesx = ndivx + 1; nodesy = ndivy+1; nodesz = ndivz+1

    # Obtain number of nodes and elements (5 tetrahedra per cube)
    n_nodes = (ndivx+1)*(ndivy+1)*(ndivz+1)
    n_elements = ndivx*ndivy*ndivz*5

    # Create nodes
    v_nodes = [zeros(3) for _ = 1:n_nodes]
    n = 1
    for i = 0:ndivx
        for j = 0:ndivy
            for k = 0:ndivz
                v_nodes[n] = [i*gapx, j*gapy, k*gapz]
                n += 1
            end
        end
    end

    # Create mesh
    v_mesh = []
    v_mesh = [zeros(Int64, 4) for _ = 1:n_elements] # This array will store the nodes pertaining to each tetrahedron. Each row contains the numbers assigned to the nodes of a tetrahedron in it. Therefore length(v_mesh) = elements
    index = 1
    for i = 1:ndivx
        for j = 1:ndivy
            for k = 1:ndivz
                ref_pt = (i-1)*nodesy*nodesz + (j-1)*nodesz + k # Reference point of each generating cube (its nearest point to the origin)

                # Each line is one of the 5 tetrahedrons within a generating cube
                # Nodal order of each tetrahedron is chosen according to OptiStruct documentation

                v_mesh[index+4] = [ref_pt + nodesy*nodesz,  ref_pt + nodesz,  ref_pt + 1,  ref_pt + nodesy*nodesz + nodesz + 1]
                v_mesh[index+1] = [ref_pt + nodesy*nodesz,  ref_pt + nodesy*nodesz + nodesz,  ref_pt + nodesz,  ref_pt + nodesy*nodesz + nodesz + 1]
                v_mesh[index+3] = [ref_pt + nodesy*nodesz,  ref_pt + nodesy*nodesz + nodesz + 1,  ref_pt + 1,  ref_pt + nodesy*nodesz + 1]
                v_mesh[index+2] = [ref_pt + 1,  ref_pt + nodesy*nodesz + nodesz + 1,  ref_pt + nodesz, ref_pt + nodesz + 1]
                v_mesh[index] = [ref_pt + nodesy*nodesz,  ref_pt + nodesz,  ref_pt,  ref_pt + 1]

                index += 5
            end
        end
    end
    # for i in range(1, stop=ndivx)
    #     for j in range(1, stop=ndivy)
    #         for k in range(1, stop=ndivz)
    #             global n_glob
    #             # numbering: i+j*nx+k*nx*ny
    #             # [i,j,k] es la posicion que ocupa el nodo local n(1,2,3...) en la malla
    #             n1=[i,j,k]
    #             n2=[i+1,j,k]
    #             n4=[i,j+1,k]
    #             n3=[i+1,j+1,k]
    #             n5=[i,j,k+1]
    #             n6=[i+1,j,k+1]
    #             n7=[i+1,j+1,k+1]
    #             n8=[i,j+1,k+1]

    #             n_loc=[n1,n2,n3,n4,n5,n6,n7,n8]
    #             n_glob=[]
    #             for l in range(1, stop=length(n_loc))
    #                 #n_glob.append(n_loc[l][0]*(ndivx+1)*(ndivy+1)+n_loc[l][1]*(ndivz+1)+n_loc[l][2])
    #                 push!(n_glob, (n_loc[l][1]-1)*(ndivz+1)*(ndivy+1)+(n_loc[l][2]-1)*(ndivz+1)+(n_loc[l][3]-1) + 1)
    #             end

    #             # El vec_aux nos va a permitir obtener los tetraedros a partir de los elementos cubicos.
    #             vec_aux=[[1,2,4,5],[2,3,4,7],[5,8,7,4],[5,7,6,2],[5,7,2,4]]
    #             for l in range(1,stop=length(vec_aux))
    #                 push!(v_mesh, [n_glob[vec_aux[l][1]], n_glob[vec_aux[l][2]], n_glob[vec_aux[l][3]], n_glob[vec_aux[l][4]]])
    #             end
    #         end
    #     end
    # end

    ### WRITE and SAVE
    # Check if mesh-data dir exists, creates it otherwise
    if !isdir("mesh-data")
        mkdir("mesh-data")
    end

    fini= open("./mesh-data/mesh.txt","w")
    fmech= open("./mesh-data/mechanical_mat.txt","w")
    mechanical_mesh = [zeros(2) for _ = 1:n_elements]

    for i in range(1,stop=length(v_mesh))
        for j in range(1,stop=length(v_mesh[i]))
            write(fini, string((v_mesh[i][j])))
            write(fini, "  ")
        end
        if i < length(v_mesh)
            write(fini, "\n")
        end

        mechanical_mesh[i] = [E, ν]
        write(fmech, string( E  ))
        write(fmech, "  ")
        write(fmech, string( ν  ))
        write(fmech, "\n")
    end

    close(fini)
    close(fmech)

    return (n_nodes, n_elements, v_nodes, v_mesh, mechanical_mesh)
end

end #module
