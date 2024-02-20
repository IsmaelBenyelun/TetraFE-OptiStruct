module OptiStruct

export optistruct

"""
    u, strains_el, stresses_el = optistruct(fileName)
Run OptiStruct to solve the FE problem. From the files generated by OptiStruct; read displacements, strains and stresses in each element
"""
function optistruct(fileName::String, solver_path::String, dof_node::Int64)

    # Run the .fem file in OptiStruct. Change directories if necessary
    fem_path = ".\\OptistructFiles\\" * fileName * ".fem" # Directory of the .fem file

    run(Cmd(`$solver_path $fem_path -optskip \> NUL`)) # Optskip analysis is performed by OptiStruct


    ###
    ## If no .strn, .strs. or .disp files are created by OptiStruct, an error must have ocurred during running. Check .out file to find the error.
    ###

    ## Read data from the output files

    # Displacements
    open(".\\OptistructFiles\\" * fileName * ".disp", "r") do disp_read
        global disp_lines = readlines(disp_read)
    end

    u = zeros(dof_node * length(disp_lines[3:end])) # First 2 lines of the file show data related to the subcase that has been calculated, not displacements
    i = 1
    for line in disp_lines[3:end] # Each line is a string that contains the displacements of a node.

        u_x = parse(Float64, strip(line[11:22]))
        u_y = parse(Float64, strip(line[23:34]))
        u_z = parse(Float64, strip(line[35:46]))

        u[i : i+2] = [u_x, u_y, u_z]
        i += 3
    end


    # Strains
    open(".\\OptistructFiles\\" * fileName * ".strn", "r") do strn_read
        global strains_lines = readlines(strn_read)
    end

    # strains_el = [zeros((3,3)) for _ = 1:length(strains_lines[3:end])] # First 2 lines of the file show data related to the subcase that has been calculated, not strains
    strains_el = [zeros(6) for _ = 1:length(strains_lines[3:end])] # First 2 lines of the file show data related to the subcase that has been calculated, not strains
    i = 1
    for line in strains_lines[3:end] # Each line is a string that contains the strains of an element

        eps_xx = parse(Float64, strip(line[23:34]))
        eps_yy = parse(Float64, strip(line[35:46]))
        eps_zz = parse(Float64, strip(line[47:58]))
        eps_xy = parse(Float64, strip(line[59:70]))
        eps_yz = parse(Float64, strip(line[71:82]))
        eps_xz = parse(Float64, strip(line[83:94]))

        # strains_el[i] = [eps_xx eps_xy eps_xz; eps_xy eps_yy eps_yz; eps_xz eps_yz eps_zz] # matrix form
        strains_el[i] = [eps_xx,  eps_yy, eps_zz, eps_xy, eps_yz, eps_xz]
        i += 1
    end


    # Stresses
    open(".\\OptistructFiles\\" * fileName * ".strs", "r") do strs_read
        global stresses_lines = readlines(strs_read)
    end

    # stresses_el = [zeros((3,3)) for _ = 1:length(stresses_lines[3:end])] # First 2 lines of the file show data related to the subcase that has been calculated, not stresses
    stresses_el = [zeros(6) for _ = 1:length(stresses_lines[3:end])] # First 2 lines of the file show data related to the subcase that has been calculated, not stresses
    i = 1
    for line in stresses_lines[3:end] # Each line is a string that contains the strains of an element

        sigma_xx = parse(Float64, strip(line[23:34]))
        sigma_yy = parse(Float64, strip(line[35:46]))
        sigma_zz = parse(Float64, strip(line[47:58]))
        sigma_xy = parse(Float64, strip(line[59:70]))
        sigma_yz = parse(Float64, strip(line[71:82]))
        sigma_xz = parse(Float64, strip(line[83:94]))

        # stresses_el[i] = [sigma_xx sigma_xy sigma_xz; sigma_xy sigma_yy sigma_yz; sigma_xz sigma_yz sigma_zz]
        stresses_el[i] = [sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_xz]
        i += 1
    end

    ip_ID = [zeros(Int64, 4) for _ = 1:length(strains_lines[3:end])] # Arbitrary value. Will never be used when OptiStruct is the solver. Defining ip_ID is necessary to use multiple dispatch in energy_distribution and H_generator functions

    return u, strains_el, stresses_el, ip_ID
end

end
