module ParaView

export PrintParaview

function PrintParaview(v_nodes::Array{Array{Float64,1}, 1}, v_mesh::Array{Array{Int64, 1},1}, name::String)
    #v_nodes es el array que contiene las coordenadas de cada nodo
    #propiedades es el array que contiene los resultados de cada nodo: stress, strains...
    #v_mesh es el array que indica qué nodos pertenecen a qué elemento
    num_pts=length(v_nodes[:,1]) # 3D. #Esto es igual al numero de nodos de la estructura

    dim_tet=length(v_mesh[:,1]) # Tetrahedro. Esto es igual al numero de elementos de la estructura

    nam = name

    # v_tetra should be a vector, defined previously
    #f=open("output_micro.vtu","w")
    f=open(nam,"w")
    # la funcion \n implica un salto de linea en el archivo de escritura
    write(f, """<?xml version="1.0"?> \n""")
    write(f, """<VTKFile type="UnstructuredGrid"> \n""")
    write(f, "\t<UnstructuredGrid> \n")
    write(f, """\t\t<Piece NumberOfPoints=" """)
    write(f, string(Int(num_pts)))
    write(f, """ "  NumberOfCells=" """)
    write(f, string(Int(dim_tet)))
    write(f, """ ">  \n""")
    write(f, "\t\t\t<Points>  \n")
    write(f, """\t\t\t\t<DataArray type="Float64" NumberOfComponents="3" format="ascii">   \n""")

    #Data array --> Coordinates
    #Para cada uno de los nodos (contador i) escribimos sus coordenadas
    for i in range(1,stop=num_pts)
        #coordenada x
        write(f, "\t\t\t\t\t", string(v_nodes[i][1]))
        write(f, "   ")
        #coordenada y
        write(f, string(v_nodes[i][2]))
        write(f, "   ")
        #coordenada z
        write(f, string(v_nodes[i][3]))
        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>   \n")
    write(f, "\t\t\t</Points>  \n")
    write(f, "\t\t\t<Cells>  \n")
    write(f, """\t\t\t\t<DataArray type="Int32" Name="connectivity" format="ascii">   \n""")

    # Data array --> Connectivity
    # Bucle para cada uno de los elementos. Se va a seleccionar cada uno de los nodos
    # del elemento i.
    for i in range(1,stop=dim_tet)
        write(f, "\t\t\t\t\t", string(v_mesh[i][1]-1))
        write(f, "   ")
        write(f, string(v_mesh[i][2]-1))
        write(f, "   ")
        write(f, string(v_mesh[i][3]-1))
        write(f, "   ")
        write(f, string(v_mesh[i][4]-1))
        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>   \n")
    write(f, """\t\t\t\t<DataArray type="Int32" Name="offsets" format="ascii">   \n""")

    # Data array --> Offset
    #
    for i in range(1,stop=dim_tet)
        aux_p=(i)*4
        write(f, "\t\t\t\t\t", string(aux_p))
        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>   \n")
    write(f, """\t\t\t\t<DataArray type="UInt8" Name="types" format="ascii">   \n""")

    for i in range(1,stop=dim_tet)
        write(f, "\t\t\t\t\t", string(10))
        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>   \n")
    write(f, "\t\t\t</Cells>  \n")

    write(f, "\t\t\t<PointData>  \n")

    write(f, "\t\t\t</PointData>  \n")
    write(f, "\t\t</Piece>  \n")
    write(f, "\t</UnstructuredGrid> \n")
    write(f, "</VTKFile> \n")
    close(f)

    return;
end

function PrintParaview(v_nodes, v_mesh, strains, stresses, displacements, name::String)
    num_pts=length(v_nodes[:,1]) # 3D. #Esto es igual al numero de nodos de la estructura
    dim_tet=length(v_mesh[:,1]) # Tetrahedro. Esto es igual al numero de elementos de la estructura

    f = open(name,"w")
    # la funcion \n implica un salto de linea en el archivo de escritura
    write(f, """<?xml version="1.0"?>\n""")
    write(f, """<VTKFile type="UnstructuredGrid">\n""")
    write(f, "\t<UnstructuredGrid>\n")
    write(f, """\t\t<Piece NumberOfPoints=" """)
    write(f, string(Int(num_pts)))
    write(f, """ "  NumberOfCells=" """)
    write(f, string(Int(dim_tet)))
    write(f, """ ">\n""")
    write(f, "\t\t\t<Points>\n")
    write(f, """\t\t\t\t<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n""")

    #Data array --> Coordinates
    #Para cada uno de los nodos (contador i) escribimos sus coordenadas
    for i in range(1,stop=num_pts)
        #coordenada x
        write(f, "\t\t\t\t\t", string(v_nodes[i][1]))
        write(f, "   ")
        #coordenada y
        write(f, string(v_nodes[i][2]))
        write(f, "   ")
        #coordenada z
        write(f, string(v_nodes[i][3]))
        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>\n")
    write(f, "\t\t\t</Points>\n")
    write(f, "\t\t\t<Cells>\n")
    write(f, """\t\t\t\t<DataArray type="Int32" Name="connectivity" format="ascii">\n""")

    # Data array --> Connectivity
    # Bucle para cada uno de los elementos. Se va a seleccionar cada uno de los nodos
    # del elemento i.
    for i in range(1,stop=dim_tet)
        write(f, "\t\t\t\t\t", string(v_mesh[i][1]-1))
        write(f, "   ")
        write(f, string(v_mesh[i][2]-1))
        write(f, "   ")
        write(f, string(v_mesh[i][3]-1))
        write(f, "   ")
        write(f, string(v_mesh[i][4]-1))
        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>\n")
    write(f, """\t\t\t\t<DataArray type="Int32" Name="offsets" format="ascii">\n""")

    # Data array --> Offset
    #
    for i in range(1,stop=dim_tet)
        aux_p=(i)*4
        write(f, "\t\t\t\t\t", string(aux_p))
        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>\n")
    write(f, """\t\t\t\t<DataArray type="UInt8" Name="types" format="ascii">\n""")

    for i in range(1,stop=dim_tet)
        write(f, "\t\t\t\t\t", string(10))
        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>\n")
    write(f, "\t\t\t</Cells>\n")

    write(f, "\t\t\t<PointData>\n")
    # Este apartado es el que se debe copiar y pegar para cada una de las propiedades
    # que queremos estudiar.

    write(f, """\t\t\t\t<DataArray type="Float64" Name="Displacements" NumberOfComponents="3" format="ascii">\n""")

    for i in range(1,stop=num_pts)
        write(f,  "\t\t\t\t\t", string(displacements[i*3-2]))
        write(f, "   ")
        write(f, string(displacements[i*3-1]))
        write(f, "   ")
        write(f, string(displacements[i*3]))

        write(f, "\n")
    end

    write(f, "\t\t\t\t</DataArray>\n")
    write(f, "\t\t\t</PointData>\n")

    write(f, "\t\t\t<CellData>\n")
    write(f, """\t\t\t\t<DataArray type="Float64" Name="Strains" NumberOfComponents="6" format="ascii">\n""")

    for i in range(1,stop=dim_tet)
        write(f, "\t\t\t\t\t", string(strains[i][1]))
        write(f, "   ")
        write(f,  string(strains[i][2]))
        write(f, "   ")
        write(f,  string(strains[i][3]))
        write(f, "   ")
        write(f,  string(strains[i][4]))
        write(f, "   ")
        write(f,  string(strains[i][5]))
        write(f, "   ")
        write(f,  string(strains[i][6]))

        write(f, "\n")
    end
    write(f, "\t\t\t\t</DataArray>\n")
    #
    write(f, """\t\t\t\t<DataArray type="Float64" Name="Stresses" NumberOfComponents="6" format="ascii">\n""")

    for i in range(1,stop=dim_tet)
        write(f, "\t\t\t\t\t", string(stresses[i][1]))
        write(f, "   ")
        write(f,  string(stresses[i][2]))
        write(f, "   ")
        write(f,  string(stresses[i][3]))
        write(f, "   ")
        write(f,  string(stresses[i][4]))
        write(f, "   ")
        write(f,  string(stresses[i][5]))
        write(f, "   ")
        write(f,  string(stresses[i][6]))

        write(f, "\n")
    end
    write(f, "\t\t\t\t</DataArray>\n")
    #

    write(f, "\t\t\t</CellData>\n")
    write(f, "\t\t</Piece>\n")
    write(f, "\t</UnstructuredGrid>\n")
    write(f, "</VTKFile>\n")
    close(f)

    return;
end

end
