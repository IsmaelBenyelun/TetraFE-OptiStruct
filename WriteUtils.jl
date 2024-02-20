module WriteUtils

using CSV
using DataFrames
using DataStructures
#
using MyStructs
export write_input_file
"""
    fileName = write_input_file(solver, material, v_nodes, v_mesh, Properties_mesh, dof_node, u, f, restricted_DOF, force_DOF)

Write the .fem/.bdf file which will be read by OptiStruct/Nastran. The .fem/.bdf file will be saved in a new folder, where all OptiStruct/Nastran generated files will be saved.
"""
function write_input_file(solver::String, elem_type::element_type, material::material_type, v_nodes::Array{Array{Float64, 1}, 1}, v_mesh::Array{Array{Int64, 1}, 1}, Properties_mesh::Array{Array{Float64, 1}, 1}, dof_node::Int64, u::Array{Float64, 1}, f::Array{Float64, 1},
    restricted_DOF::Array{Int64, 1}, force_DOF::Array{Int64, 1})

    # Create a list of dictionaries with the nodal data
    GRID_data = grid_data(v_nodes)

    # Create a list of dictionaries with the element data
    ELEMENT_data = element_data(elem_type, v_mesh)

    # Create a list of dictionaries with the material data of each element
    MATERIAL_data = material_data(solver, material, Properties_mesh)

    # Create a list of dictionaries with the porperty data of each element
    PROPERTY_data = property_data(solver, Properties_mesh)

    # Create a list of dictionaries with the boundary conditions and restrictions data
    SPC_data = boundary_conditions_data(dof_node, u, restricted_DOF)

    # Create a list of dictionaries with the force data
    FORCE_data = force_data(dof_node, f, force_DOF)

    # Create the folder where the files are going to be saved. If already existent, remove its files
    folderName = titlecase(solver) * "Files"
    if isdir(folderName)
        foreach(rm, readdir(folderName, join=true))
    else
        mkdir(folderName)
    end

    # Output to a .fem/.bdf file. File is generated in the new folder
    if solver == "OPTISTRUCT"

        fileName = write_data_optistruct(elem_type, material, GRID_data, ELEMENT_data, MATERIAL_data, PROPERTY_data, SPC_data, FORCE_data)

    elseif solver == "NASTRAN"

        fileName = write_data_nastran(elem_type, material, GRID_data, ELEMENT_data, MATERIAL_data, PROPERTY_data, SPC_data, FORCE_data)

    end

    return fileName
end



"""
    GRID_data = grid_data(v_nodes)

Store the coordinates of the nodes in dictionaries. A dictionary is created for each node. Therefore, GRID_data is a list of (no. nodes) dictionaries
"""
function grid_data(v_nodes::Array{Array{Float64, 1}, 1})

    GRID_data = Vector{OrderedDict}()
    ID = 1

    for node in v_nodes
        push!(GRID_data, OrderedDict(
            [("label", "GRID"), ("ID", ID), ("CP", ""), ("X1", node[1]), ("X2", node[2]), ("X3", node[3]), ("CD", ""), ("PS", "")]
            )
        )
        ID += 1
    end

    return GRID_data
end



"""
    ELEMENT_data = element_data(elem_type, v_mesh)

Tetra4 element
Store the nodes which belong to each element in dictionaries. A dictionary is created for each element. Therefore, ELEMENT_data is a list of (no. elements) dictionaries
"""
function element_data(elem_type::Tetra4, v_mesh::Array{Array{Int64, 1}, 1})

    ELEMENT_data = Vector{OrderedDict}()
    ID = 1

    for element in v_mesh
        push!(ELEMENT_data, OrderedDict(
            [("label", "CTETRA"), ("EID", ID), ("PID", ID), ("G1", element[1]), ("G2", element[2]), ("G3", element[3]), ("G4", element[4])]
            )
        )
        ID += 1
    end

    return ELEMENT_data
end


"""
    (ELEMENT_data1, ELEMENT_data2) = element_data(elem_type, v_mesh)

Tetra10 element
Store the nodes which belong to each element in dictionaries. A dictionary is created for each element. Therefore, ELEMENT_data is a list of (no. elements) dictionaries
"""
function element_data(elem_type::Tetra10, v_mesh::Array{Array{Int64, 1}, 1})

    ELEMENT_data1 = Vector{OrderedDict}(); ELEMENT_data2 = Vector{OrderedDict}()
    ID = 1

    for element in v_mesh
        push!(ELEMENT_data1, OrderedDict(
            [("label", "CTETRA"), ("EID", ID), ("PID", ID), ("G1", element[1]), ("G2", element[2]), ("G3", element[3]), ("G4", element[4]), ("G5", element[5]), ("G6", element[6])]
            )
        )
        push!(ELEMENT_data2, OrderedDict([("extraln", "+"), ("G7", element[7]), ("G8", element[8]), ("G9", element[9]), ("G10", element[10])]))
        ID += 1
    end

    return ELEMENT_data1, ELEMENT_data2
end



"""
    MATERIAL_data = material_data(solver, material, Properties_mesh)

Isotropic material
Store mechanical properties of the materials in dictionaries. A dictionary is created for each material. Therefore, since each element has different material, MATERIAL_data is a list of (no. elements) dictionaries
"""
function material_data(solver::String, material::Isotropic, Properties_mesh::Array{Array{Float64, 1}, 1})

    MATERIAL_data = Vector{OrderedDict}()
    ID = 1

    for property in Properties_mesh
        push!(MATERIAL_data, OrderedDict([("label", "MAT1"), ("MID", ID), ("E", property[1]), ("G", ""), ("NU", property[2])]))
        ID += 1
    end

    return MATERIAL_data
end



"""
    (MATERIAL_data1, MATERIAL_data2) = material_data(solver, material, Properties_mesh)

Orthotropic material
Material data is stored in dictionaries
To write using CSV.write afterwards, is more convenient to store the mechanical properties of each element in 2 different lists of (no.elements)
MATERIAL_data1 stores the first line of each material data (E1, E2, E3, nu12, nu23, nu31) and MATERIAL_data2 stores the second line (G12, G23, G31)
"""
function material_data(solver::String, material::Orthotropic, Properties_mesh::Array{Array{Float64, 1}, 1})

    if solver == "OPTISTRUCT"

        label = "MAT9OR"

    elseif solver == "NASTRAN"

        label = "MATORT"
    end

    MATERIAL_data1 = Vector{OrderedDict}()
    MATERIAL_data2 = Vector{OrderedDict}()
    ID = 1

    for property in Properties_mesh

        # Round properties to avoid line lengths greater than 80 characters in .fem/.bdf file
        E1 = round(property[1]; digits=4)
        E2 = round(property[2]; digits=4)
        E3 = round(property[3]; digits=4)
        NU12 = round(property[4]; digits=7)
        NU23 = round(property[5]; digits=7)
        NU31 = round(property[3]/property[1] * property[6]; digits=7)

        push!(MATERIAL_data1, OrderedDict(
            [("label", label), ("MID", ID), ("E1", E1), ("E2", E2), ("E3", E3), ("NU12", NU12), ("NU23", NU23), ("NU31", NU31)]
            )
        )
        push!(MATERIAL_data2, OrderedDict([("extraln", "+"), ("G12", property[7]), ("G23", property[8]), ("G31", property[9])]))
        ID += 1
    end

    return MATERIAL_data1, MATERIAL_data2
end



"""
    PROPERTY_data = property_data(solver, Properties_mesh)

Property data is stored in dictionaries. A dictionary is created for each element because each element has different material
"""
function property_data(solver::String, Properties_mesh::Array{Array{Float64, 1}, 1})

    PROPERTY_data = Vector{OrderedDict}()
    ID = 1

    if solver == "OPTISTRUCT"

        for ID = 1 : length(Properties_mesh)
            push!(PROPERTY_data, OrderedDict([("label", "PSOLID"), ("PID", ID), ("MID1", ID), ("CORDM", "")]))
            ID += 1
        end

    elseif solver == "NASTRAN"

        for ID = 1 : length(Properties_mesh)
            push!(PROPERTY_data, OrderedDict([("label", "PSOLID"), ("PID", ID), ("MID", ID), ("CORDM", ""), ("IN", ""), ("STRESS", "GAUSS"), ("ISOP", "")]))
            ID += 1
        end

    end

    return PROPERTY_data
end



"""
    SPC_data = boundary_conditions(dof_node, u, restricted_DOF)

Boundary conditions and imposed displacements. A dictionary is created for each SPC. Therefore, SPC_data is a list of (no. SPC) dictionaries
"""
function boundary_conditions_data(dof_node::Int64, u::Array{Float64, 1}, restricted_DOF::Array{Int64, 1})

    SPC_data = Vector{OrderedDict}()
    ID = 1

    for DOF in restricted_DOF

        GID = (DOF-1) ÷ dof_node + 1

        nodal_DOF = DOF - dof_node*GID + dof_node

        push!(SPC_data, OrderedDict([("label", "SPC"), ("SID", ID), ("GID", GID), ("C", nodal_DOF), ("D", u[DOF])]))
    end

    return SPC_data
end



"""
    FORCE_data = force_data(dof_node, f, force_DOF)

Force. A dictionary is created for each node force. Therefore, FORCE_data is a list of (no. forces) dictionaries
"""
function force_data(dof_node::Int64, f::Array{Float64, 1}, force_DOF::Array{Int64, 1})

    FORCE_data = Vector{OrderedDict}()
    ID = 1

    for DOF in force_DOF

        GID = (DOF-1) ÷ dof_node + 1

        if DOF % dof_node == 1
            push!(FORCE_data, OrderedDict([("label", "FORCE"), ("SID", ID), ("GID", GID), ("CID", ""), ("F", f[DOF]), ("N1", 1.), ("N2", 0.), ("N3", 0.)]))

        elseif DOF % dof_node == 2
            push!(FORCE_data, OrderedDict([("label", "FORCE"), ("SID", ID), ("GID", GID), ("CID", ""), ("F", f[DOF]), ("N1", 0.), ("N2", 1.), ("N3", 0.)]))

        elseif DOF % dof_node == 0
            push!(FORCE_data, OrderedDict([("label", "FORCE"), ("SID", ID), ("GID", GID), ("CID", ""), ("F", f[DOF]), ("N1", 0.), ("N2", 0.), ("N3", 1.)]))

        end
    end

    return FORCE_data
end



"""
    fileName = write_data(material, GRID_data, ELEMENT_data, MATERIAL_data, PORPERTY_data, SPC_data, FORCE_data)

Isotropic case
Write all the data in a .fem file, following the order prescribed by OptiStruct/Nastran
"""
function write_data_optistruct(elem_type::element_type, material::Isotropic, GRID_data::Array{OrderedDict, 1}, ELEMENT_data::Union{Array{OrderedDict, 1}, Tuple{Array{OrderedDict, 1}, Array{OrderedDict, 1}}},
    MATERIAL_data::Array{OrderedDict, 1}, PROPERTY_data::Array{OrderedDict, 1}, SPC_data::Array{OrderedDict, 1}, FORCE_data::Array{OrderedDict, 1})

    fileName = "3d_tetra"

    open(".\\OptistructFiles\\" * fileName * ".fem", "w") do fem

        # I/O Options Section
        write(fem, "OUTPUT, H3D, NONE\n") # Prevent OptiStruct from generating .h3d file
        write(fem, "OUTPUT, HM, NONE\n") # Prevent OptiStruct from generating .hm file
        write(fem, "OUTPUT, MVW, NO\n") # Prevent OptiStruct from generating .mvw file

        # Subcase Information Section
        write(fem, "SUBCASE,1\n")
        write(fem, "  LABEL loadstep1\n")
        write(fem, "ANALYSIS STATICS\n") # Perform linear static analysis
        write(fem, "  SPC = 1\n")
        write(fem, "  LOAD = 1\n")
        write(fem, "  DISPLACEMENT(OPTI) = ALL\n") # Request displacement vector output in OptiStruct results format (.disp)
        write(fem, "  STRAIN(OPTI) = ALL\n") # Request strain output in OptiStruct results format (.strns)
        write(fem, "  STRESS(OPTI) = ALL\n\n") # Request stress output in OptiStruct results format (.strs)

        # Bulk Data Section
        write(fem, "BEGIN BULK\n\n")

        CSV.write(fem, DataFrame(GRID_data), append=true); write(fem, "\n")

        if typeof(elem_type) == Tetra4

            CSV.write(fem, DataFrame(ELEMENT_data), append=true); write(fem, "\$\n")

        elseif typeof(elem_type) == Tetra10

            ELEMENT_data1 = ELEMENT_data[1]
            ELEMENT_data2 = ELEMENT_data[2]

            for i = 1 : length(ELEMENT_data1)
                CSV.write(fem, DataFrame(ELEMENT_data1[i]), append=true)
                CSV.write(fem, DataFrame(ELEMENT_data2[i]), append=true)
            end
            write(fem, "\n")
        end

        CSV.write(fem, DataFrame(PROPERTY_data), append=true); write(fem, "\n")
        CSV.write(fem, DataFrame(MATERIAL_data), append=true); write(fem, "\$\n")

        CSV.write(fem, DataFrame(SPC_data), append=true); write(fem, "\n")
        CSV.write(fem, DataFrame(FORCE_data), append=true); write(fem, "\n")

        write(fem, "ENDDATA")
    end

    return fileName
end



"""
    fileName = write_data(material, GRID_data, ELEMENT_data, MATERIAL_data1, MATERIAL_data2, PORPERTY_data, SPC_data, FORCE_data)

Orthotropic case
Write all the data in a .fem file, following the order prescribed by OptiStruct/Nastran
"""
function write_data_optistruct(elem_type::element_type, material::Orthotropic, GRID_data::Array{OrderedDict, 1}, ELEMENT_data::Union{Array{OrderedDict, 1}, Tuple{Array{OrderedDict, 1}, Array{OrderedDict, 1}}},
    MATERIAL_data::Tuple{Array{OrderedDict, 1}, Array{OrderedDict, 1}}, PROPERTY_data::Array{OrderedDict, 1}, SPC_data::Array{OrderedDict, 1}, FORCE_data::Array{OrderedDict, 1})

    MATERIAL_data1 = MATERIAL_data[1]
    MATERIAL_data2 = MATERIAL_data[2]

    fileName = "3d_tetra"

    open(".\\OptistructFiles\\" * fileName * ".fem", "w") do fem

        # I/O Options Section
        write(fem, "OUTPUT, H3D, NONE\n") # Prevent OptiStruct from generating .h3d file
        write(fem, "OUTPUT, HM, NONE\n") # Prevent OptiStruct from generating .hm file
        write(fem, "OUTPUT, MVW, NO\n") # Prevent OptiStruct from generating .mvw file

        # Subcase Information Section
        write(fem, "SUBCASE,1\n")
        write(fem, "  LABEL loadstep1\n")
        write(fem, "ANALYSIS STATICS\n") # Perform linear static analysis
        write(fem, "  SPC = 1\n")
        write(fem, "  LOAD = 1\n")
        write(fem, "  DISPLACEMENT(OPTI) = ALL\n") # Request displacement vector output in OptiStruct results format
        write(fem, "  STRAIN(OPTI) = ALL\n") # Request strain output in OptiStruct results format
        write(fem, "  STRESS(OPTI) = ALL\n\n") # Request stress output in OptiStruct results format

        # Bulk Data Section
        write(fem, "BEGIN BULK\n\n")

        CSV.write(fem, DataFrame(GRID_data), append=true); write(fem, "\n")

        if typeof(elem_type) == Tetra4
            CSV.write(fem, DataFrame(ELEMENT_data), append=true); write(fem, "\$\n")

        elseif typeof(elem_type) == Tetra10

            ELEMENT_data1 = ELEMENT_data[1]
            ELEMENT_data2 = ELEMENT_data[2]

            for i = 1 : length(ELEMENT_data1)
                CSV.write(fem, DataFrame(ELEMENT_data1[i]), append=true)
                CSV.write(fem, DataFrame(ELEMENT_data2[i]), append=true)
            end
            write(fem, "\n")
        end

        CSV.write(fem, DataFrame(PROPERTY_data), append=true); write(fem, "\n")

        for i = 1 : length(MATERIAL_data1)
            CSV.write(fem, DataFrame(MATERIAL_data1[i]), append=true)
            CSV.write(fem, DataFrame(MATERIAL_data2[i]), append=true)
        end
        write(fem, "\n")

        CSV.write(fem, DataFrame(SPC_data), append=true); write(fem, "\n")
        CSV.write(fem, DataFrame(FORCE_data), append=true); write(fem, "\n")

        write(fem, "ENDDATA")
    end

    return fileName
end



"""
    fileName = write_data_nastran(material, GRID_data, ELEMENT_data, MATERIAL_data, PORPERTY_data, SPC_data, FORCE_data)

Isotropic case
Write all the data in a .bdf file, following the order prescribed by Nastran
"""
function write_data_nastran(elem_type::element_type, material::Isotropic, GRID_data::Array{OrderedDict, 1}, ELEMENT_data::Union{Array{OrderedDict, 1}, Tuple{Array{OrderedDict, 1}, Array{OrderedDict, 1}}},
    MATERIAL_data::Array{OrderedDict, 1}, PROPERTY_data::Array{OrderedDict, 1}, SPC_data::Array{OrderedDict, 1}, FORCE_data::Array{OrderedDict, 1})

    fileName = "3d_tetra"

    open(".\\NastranFiles\\" * fileName * ".bdf", "w") do bdf

        # Executive Control Statements
        write(bdf, "SOL 400\n") # Execute SOL 400 sequence
        write(bdf, "CEND\n")

        write(bdf, "\$\n")

        # Case Control Commands
        write(bdf, "ECHO = NONE\n") # Neither sorted nor unsorted Bulk Data will be printed
        write(bdf, "SUBCASE = 1\n")
        write(bdf, "  STEP = 1\n")
        write(bdf, "    ANALYSIS = LNSTATICS\n") # Perform linear static analysis
        write(bdf, "    SPC = 2\n")
        write(bdf, "    LOAD = 2\n")
        write(bdf, "    DISPLACEMENT(PUNCH) = ALL\n") # Request displacement vector in a punch file (.pch)
        write(bdf, "    STRAIN(PUNCH) = ALL\n") # Request strain output in a punch file (.pch)
        write(bdf, "    STRESS(PUNCH) = ALL\n") # Request stress output in a punch file (.pch)
        write(bdf, "BEGIN BULK\n")

        write(bdf, "\$\n")

        # Bulk Data Entries
        write(bdf, "PARAM,PRTMAXIM,YES\n") # Printout  maximums of applied loads, single-point forces of constraint, multipoint forces of constraint, and displacements

        write(bdf, "\$\n")

        CSV.write(bdf, DataFrame(GRID_data), append=true); write(bdf, "\$\n")

        if typeof(elem_type) == Tetra4
            CSV.write(bdf, DataFrame(ELEMENT_data), append=true); write(bdf, "\$\n")

        elseif typeof(elem_type) == Tetra10

            ELEMENT_data1 = ELEMENT_data[1]
            ELEMENT_data2 = ELEMENT_data[2]

            for i = 1 : length(ELEMENT_data1)
                CSV.write(bdf, DataFrame(ELEMENT_data1[i]), append=true)
                CSV.write(bdf, DataFrame(ELEMENT_data2[i]), append=true)
            end
            write(bdf, "\$\n")
        end

        CSV.write(bdf, DataFrame(PROPERTY_data), append=true); write(bdf, "\$\n")
        CSV.write(bdf, DataFrame(MATERIAL_data), append=true); write(bdf, "\$\n")

        write(bdf, "SPCADD,2,1\n") # Define a single-point constraint set as a union of single-point constraint sets defined on SPC
        write(bdf, "LOAD,2,1.,1.,1\n") # Defines a static load as a linear combination of load sets defined via FORCE
        write(bdf, "\$\n")
        CSV.write(bdf, DataFrame(SPC_data), append=true); write(bdf, "\$\n")
        CSV.write(bdf, DataFrame(FORCE_data), append=true); write(bdf, "\$\n")

        write(bdf, "ENDDATA")

    end

    return fileName
end



"""
    fileName = write_data(material, GRID_data, ELEMENT_data, MATERIAL_data1, MATERIAL_data2, PORPERTY_data, SPC_data, FORCE_data)

Orthotropic case
Write all the data in a .bdf file, following the order prescribed by Nastran
"""
function write_data_nastran(elem_type::element_type, material::Orthotropic, GRID_data::Array{OrderedDict, 1}, ELEMENT_data::Union{Array{OrderedDict, 1}, Tuple{Array{OrderedDict, 1}, Array{OrderedDict, 1}}},
    MATERIAL_data::Tuple{Array{OrderedDict, 1}, Array{OrderedDict, 1}}, PROPERTY_data::Array{OrderedDict, 1}, SPC_data::Array{OrderedDict, 1}, FORCE_data::Array{OrderedDict, 1})

    MATERIAL_data1 = MATERIAL_data[1]
    MATERIAL_data2 = MATERIAL_data[2]

    fileName = "3d_tetra"

    open(".\\NastranFiles\\" * fileName * ".bdf", "w") do bdf

        # Executive Control Statements
        write(bdf, "SOL 400\n") # Execute SOL 400 sequence
        write(bdf, "CEND\n")

        write(bdf, "\$\n")

        # Case Control Commands
        write(bdf, "ECHO = NONE\n") # Neither sorted nor unsorted Bulk Data will be printed
        write(bdf, "SUBCASE = 1\n")
        write(bdf, "  STEP = 1\n")
        write(bdf, "    ANALYSIS = LNSTATICS\n") # Perform linear static analysis
        write(bdf, "    SPC = 2\n")
        write(bdf, "    LOAD = 2\n")
        write(bdf, "    DISPLACEMENT(PUNCH) = ALL\n") # Request displacement vector in a punch file (.pch)
        write(bdf, "    STRAIN(PUNCH) = ALL\n") # Request strain output in a punch file (.pch)
        write(bdf, "    STRESS(PUNCH) = ALL\n") # Request stress output in a punch file (.pch)
        write(bdf, "BEGIN BULK\n")

        write(bdf, "\$\n")

        # Bulk Data Entries
        write(bdf, "MDLPRM,HDF5,0\n") # Create NH5RDB database
        write(bdf, "PARAM,PRTMAXIM,YES\n") # Printout  maximums of applied loads, single-point forces of constraint, multipoint forces of constraint, and displacements

        write(bdf, "\$\n")

        CSV.write(bdf, DataFrame(GRID_data), append=true); write(bdf, "\$\n")

        if typeof(elem_type) == Tetra4
            CSV.write(bdf, DataFrame(ELEMENT_data), append=true); write(bdf, "\$\n")

        elseif typeof(elem_type) == Tetra10

            ELEMENT_data1 = ELEMENT_data[1]
            ELEMENT_data2 = ELEMENT_data[2]

            for i = 1 : length(ELEMENT_data1)
                CSV.write(bdf, DataFrame(ELEMENT_data1[i]), append=true)
                CSV.write(bdf, DataFrame(ELEMENT_data2[i]), append=true)
            end
            write(bdf, "\$\n")
        end

        CSV.write(bdf, DataFrame(PROPERTY_data), append=true); write(bdf, "\$\n")

        for i = 1 : length(MATERIAL_data1)
            CSV.write(bdf, DataFrame(MATERIAL_data1[i]), append=true)
            CSV.write(bdf, DataFrame(MATERIAL_data2[i]), append=true)
        end

        write(bdf, "\$\n")

        write(bdf, "SPCADD,2,1\n") # Define a single-point constraint set as a union of single-point constraint sets defined on SPC
        write(bdf, "LOAD,2,1.,1.,1\n") # Defines a static load as a linear combination of load sets defined via FORCE
        write(bdf, "\$\n")
        CSV.write(bdf, DataFrame(SPC_data), append=true); write(bdf, "\$\n")
        CSV.write(bdf, DataFrame(FORCE_data), append=true); write(bdf, "\$\n")

        write(bdf, "ENDDATA")

    end

    return fileName
end
end
