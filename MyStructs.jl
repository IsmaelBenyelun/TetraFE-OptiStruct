module MyStructs

export material_type
export Isotropic
export Orthotropic

export element_type
export Tetra4
export Tetra10


abstract type material_type end # I have bookmarked every instance where these structs are added

mutable struct Isotropic <: material_type

    E :: Float64
    ν :: Float64 # This symbol is \nu, not a regular v

end

mutable struct Orthotropic <: material_type

    E1 :: Float64
    E2 :: Float64
    E3 :: Float64

    ν12 :: Float64
    ν23 :: Float64
    ν13 :: Float64

    G12 :: Float64
    G23 :: Float64
    G13 :: Float64

end


abstract type element_type end

struct Tetra4 <: element_type
end

struct Tetra10 <: element_type
end


end # module MyStructs
