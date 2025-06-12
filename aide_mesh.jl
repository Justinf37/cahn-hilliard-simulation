
# --------------------------------------------
#
#
# HELPER ROUTINES
#
#
# --------------------------------------------
function add_extrudedsphere(xc, yc, zc, r1, r2)

    tag1 = gmsh.model.occ.add_sphere(xc, yc, zc, r1)
    tag2 = gmsh.model.occ.add_sphere(xc, yc, zc, r2)

    out, _= gmsh.model.occ.cut([(3, tag1)],[(3, tag2)], 3)

    tag = out[1][2]
    tag
end


function add_divided_sphere(xc, yc, zc, r, n1)

    tags = collect(1: Int(n1))
    θ₁ = 2. * π / n1
    gmsh.model.occ.add_sphere(xc, yc, zc, r, tags[1], -π/2, π/2, θ₁)
    dimTags = (3, tags[1])
    for i in tags[2:end]
        gmsh.model.occ.rotate(dimTags, 0., 0., 0., 0., 0., 1., (i - 1) * θ₁)
    end

    objectDimTags = [(3, tags[1])]
    toolDimTags = [(3, tag) for tag in tags[2:end]]
    out, _ = gmsh.model.occ.fuse(objectDimTags, toolDimTags, 3)
    out
end


const _revolutions = ["ellipse", "sphere"]
const _cutspheres = ["semisphere", "truncatedsphere"]
"""
    struct Mesh{ T <: Number}

Stores the data needed in the creation of a gmsh nonuniform
mesh in Cartesian coordinates.

# Fields

- `name::String`: for rectangular objects either `rectangle` or `box`,
   whether in 2D or 3D, respectively. All in all, the following options
   are supported: `$(_allowed_names)``
- `size::NTuple{N, Number} where {N}`: the lengths of the mesh sides.

		For name == "rectangle"
			size = (a, b)

		For name == "box"
			size = (a, b, c)
		
		In the above two cases, a, b, c are the respective dimensions
		of the rectangle or the box.
		
		For name == "cylinder"
			size = (x, y, z, r)

		Where the first three coordinates define the axis vector of
		the cylinder and r is its radius.

		For name == "sphere"
			size = (r,)
		
		For name == "semisphere"
			size = (r,)
		
		For name == "truncatedsphere"
			size = (r, θ)
		
		Where r is the (semi)sphere radius and θ the polar angle from the
		interval [0, π].

        For name == "extrudedsphere"
            size = (r1, r2)
        
        Where r1 is the outer shell radius and r2 is the inner shell radius.

		For name == "ellipse"
			size = (rx, ry, rz)
		Where rx, ry and rz denote the principal axes of the ellipse.

        For name == "dividedsphere"
            size = (r, n1, n2)
        
        Where r is the sphere radius, n1 defines the partition in the ϕ
        angle and n2 the partition in the θ angle.

   
- `x0::NTuple{3, Number}`: the coordinates of the mesh origin, default to
  `(0., 0., 0.)`
- `dimension::Integer`: dimensionality of the problem, works for 2D and 3D and
   is automatically determined for a given name. As of now, only
	rectangles are two-dimensional and the rest are three-dimensional.
- `_mesh_function::Any`: an internal attribute, setting which method for
   mesh creation should be used depending on the dimensionality of the problem.

# Constructors

- `Mesh(name::String, size::NTuple{N, Number}, 
	x0::NTuple{3, Number} = (0.0, 0.0, 0.0))`
"""
mutable struct Mesh{T<:Number}

    name::String
    size::NTuple{N,T} where {N}
    x0::NTuple{3,T}
    dimension::Integer
    _mesh_function::Function
    savename::String

    function Mesh(name::String, size::NTuple{N,T},
        x0::NTuple{3,T}=(0.0, 0.0, 0.0)) where {N, T <: Number} 

        @assert name in _allowed_names "$(name) not yet implemented!"

        if name == "rectangle"
            _mesh_function = gmsh.model.occ.add_rectangle
            dim = 2
        elseif name == "box"
            _mesh_function = gmsh.model.occ.add_box
            dim = 3
        elseif name == "sphere"
            _mesh_function = gmsh.model.occ.add_sphere
            dim = 3
        elseif name == "extrudedsphere"

            _mesh_function = add_extrudedsphere
            dim = 3
        elseif name == "dividedsphere"

            _mesh_function = add_divided_sphere
            dim = 3

        elseif name == "semisphere"
            add_semisphere = (xc, yc, zc, r,) -> gmsh.model.occ.add_sphere(
                xc, yc, zc, r, -1, 0.0, pi / 2.0)
            _mesh_function = add_semisphere
            dim = 3

        elseif name == "truncatedsphere"

            add_truncatedsphere = (xc, yc, zc, r, θ) -> gmsh.model.occ.add_sphere(
                xc, yc, zc, r, -1, θ, π / 2.)
            _mesh_function = add_truncatedsphere
            dim = 3

        elseif name == "ellipse"

            function add_ellipse(xc, yc, zc, a, b, c)
                
                tag = gmsh.model.occ.add_sphere(xc, yc, zc, 1.)

                gmsh.model.occ.dilate((3, tag), xc, yc, zc, a, b, c)
                tag
            end

            _mesh_function = add_ellipse
            dim = 3

        elseif name == "cylinder"

            _mesh_function = gmsh.model.occ.add_cylinder
            dim = 3

        end
        new{T}(name, size, x0, dim, _mesh_function, "")
    end
end


function _popup(popup::Bool; close::Bool=false, wait::Float64=0.)

    if popup
        gmsh.fltk.initialize()
        gmsh.fltk.run()

        if close
            gmsh.fltk.wait(wait)
            gmsh.fltk.finalize()
        end
    end

end

function _write(write::Bool, path::String, filename::String)

    # create a path if it does not exist yet
    savename = joinpath(path, filename)
    if !isdir(path)
        mkpath(path)
    end
    # write the mesh to disk
    if write
        gmsh.write(savename)
    end
    return savename
end

function _format_filename(mesh::Mesh, lc::Any; optstring::String="")

    param_string = join(string.(mesh.size), "_")
    filename = "$(mesh.name)_shape_$(param_string)_lc_$(lc)_$(optstring).msh"

    return filename
end


"""

    genmesh(mesh::Mesh; lc::Number=5e-01, lcs::Number=5e-01, popup::Bool=false, 
        close::Bool=false, wait::Float64=0., 
        write::Bool=false, savedir::String=".", proglc::Bool=false, 
        distmin::Number=0.0, distmax::Number=1.0, lcmin::Number=5e-01, lcmax::Number=5e-01)

Generate a mesh and save it to disk (optional). The code
scans the destination folder; if a file to be generated already
exists inside the folder, skip the mesh generation step and prompt
a message instead.

## Arguments
- `mesh::Mesh`: The mesh object for which to generate the mesh.
- `lc::Number`: The characteristic length for the mesh elements.
- `lcs::Number`: The surface mesh size.
- `popup::Bool`: Whether to display the mesh in the gmsh GUI.
- `close::Bool`: Whether to close the gmsh GUI after generating the mesh.
- `wait::Float64`: The time to wait before closing the gmsh GUI.
- `write::Bool`: Whether to write the mesh to a file.
- `savedir::String`: The directory in which to save the mesh file.
- `proglc::Bool`: Whether to have finer meshing on the surface.
The following arguments are only relevant if `proglc==true`:
- `distmin::Number`: The distance at which to apply surface mesh size. Adjust according to
the size of the object and your needs.
- `distmax::Number`: The distance at which to apply volume mesh size. Adjust according to
the size of the object and your needs.
- `lcmin::Number`: The minimum mesh element size. Adjust according to your needs. By default,
this will control the minimum mesh size on the surface.
- `lcmax::Number`: The maximum mesh element size. Adjust according to your needs. By default,
this will default the maximum mesh size, hence the size in the volume.

## Returns
- `savename::String`: The name of the saved mesh file.

"""
function genmesh(mesh::Mesh; lc::Number=5e-01, 
    lcs::Number=5e-01, # surface mesh size
    popup::Bool=false,
    close::Bool=false,
    wait::Float64=0.,
    write::Bool=false,
    savedir::String=".",
    proglc::Bool=false, # whether to have finer meshing on surface
    distmin::Number=0.0,    # at which distance to apply surface lc
    distmax::Number=1.0,    # at which distance to apply volume lc
    lcmin::Number=5e-01,    # minimum lc
    lcmax::Number=5e-01)    # maximum lc

    

    _suffix = _set_prog_suffix(proglc, lcmin, lcmax, distmin, distmax)
    filename = _format_filename(mesh, lc; optstring=_suffix)
    # filename = "$(_filename)$(_suffix)"
    savename = _write(false, savedir, filename)

    if !isfile(joinpath(savedir, filename))

        println("Generating mesh!")
        # initialize the gmsh machinery
        gmsh.initialize()
        # add the model we are considering
        gmsh.model.add(mesh.name)
        # set the parameters of the grid
        body = mesh._mesh_function(mesh.x0..., mesh.size...)

        gmsh.model.occ.synchronize()
        center = gmsh.model.occ.addPoint(mesh.x0...)
        gmsh.model.occ.synchronize()
        gmsh.model.get_boundary((mesh.dimension, body),
            true, true, false)

        if mesh.name in _revolutions
            gmsh.model.mesh.embed(0, [center], mesh.dimension, 1)
            gmsh.model.addPhysicalGroup(0, [3], 1, "Center")
        end
        gmsh.model.occ.synchronize()
        # set the physical groups in the mesh
        # for easier representation
        _setgroups!(mesh, gmsh.model, body)

        if !proglc
            gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc)
        else
            _set_characteristic_lengths!(mesh, gmsh.model, body, lcmin, 
            lcmax, distmin, distmax)
        end
        # gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc)
        gmsh.model.mesh.generate(mesh.dimension)
        gmsh.model.mesh.optimize()

        # if popup == true, invoke the gmsh GUI
        savename = _write(write, savedir, filename)
        gmsh.model.mesh.getElements()
        _popup(popup, wait=wait, close=close)
        gmsh.finalize()
        println("Mesh generation done!")
    else
        println("File $(filename) already exists in folder $(savedir)")
    end

    mesh.savename = savename

end



"""
Internal function for setting physical groups of the mesh
depending on the chosen dimension of the grid.
"""
function _setgroups!(mesh::Mesh, model::Any, body::Any)

    dim = mesh.dimension
    name = mesh.name
    if name == "rectangle"

        model.addPhysicalGroup(dim - 1, [4], 1, "Left edge")
        model.addPhysicalGroup(dim - 1, [2], 2, "Right edge")
        model.addPhysicalGroup(dim - 1, [1, 2], 3, "Casing")
        model.addPhysicalGroup(dim, [body], 1, "Surface")

    elseif name == "box"

        model.addPhysicalGroup(dim - 1, [-1], 1, "Left edge")
        model.addPhysicalGroup(dim - 1, [2], 2, "Right edge")
        model.addPhysicalGroup(dim - 1, [-3, 4, -5, 6], 3, "Casing")
        model.addPhysicalGroup(dim, [body], 1, "Volume")

    elseif name in _revolutions # sphere, ellipse

        model.addPhysicalGroup(dim - 1, [-1], 1, "Casing")

        model.addPhysicalGroup(dim - 2, [1, 2, 3], 1, "Casing")
        model.addPhysicalGroup(dim- 3, [1, 2], 3, "Casing")
        model.addPhysicalGroup(dim, [body], 1, "Volume")

    elseif name == "extrudedsphere" # sphere, ellipse

        model.addPhysicalGroup(dim - 1, [1], 1, "Casing")
        model.addPhysicalGroup(dim - 2, [1, 2, 3], 2, "Casing")
        model.addPhysicalGroup(dim- 3, [1, 2], 3, "Casing")

        model.addPhysicalGroup(dim - 1, [2], 4, "Center")
        model.addPhysicalGroup(dim - 2, [4, 5, 6], 5, "Center")
        model.addPhysicalGroup(dim- 3, [3, 4], 6, "Center")
        
        model.addPhysicalGroup(dim, [body], 1, "Volume")

    elseif name in _cutspheres # semisphere, truncatedsphere
        model.addPhysicalGroup(dim - 1, [1], 1, "Casing")
        model.addPhysicalGroup(dim - 1, [2], 2, "Bottom edge")
        model.addPhysicalGroup(dim, [body], 1, "Volume")


    elseif name == "cylinder"

        model.addPhysicalGroup(dim - 1, [3], 1, "Left edge")
        model.addPhysicalGroup(dim - 1, [2], 2, "Right edge")
        model.addPhysicalGroup(dim - 1, [1], 3, "Casing")
        model.addPhysicalGroup(dim, [body], 1, "Volume")

    end

    model.occ.synchronize()

end


function _set_characteristic_lengths!(mesh::Mesh, model::Any, body::Any, lcmin, lcmax,
    distmin, distmax)

    # this is intended for 3d simulations for now
    dim = mesh.dimension
    # get entities for casing (1 is the tag of the casing)
    surface_tags = model.getEntitiesForPhysicalGroup(dim-1, 1)

    # adjust characteristic lengths on surface and elsewhere progressively
    distance_field = model.mesh.field.add("Distance")
    model.mesh.field.setNumbers(distance_field, "SurfacesList", surface_tags)

    threshold_field = gmsh.model.mesh.field.add("Threshold")
    model.mesh.field.setNumber(threshold_field, "IField", distance_field)
    model.mesh.field.setNumber(threshold_field, "LcMin", lcmin)  # Minimum element size
    model.mesh.field.setNumber(threshold_field, "LcMax", lcmax)  # Maximum element size
    model.mesh.field.setNumber(threshold_field, "DistMin", distmin )  # Distance at which SizeMin is applied
    model.mesh.field.setNumber(threshold_field, "DistMax", distmax)  # Distance at which SizeMax is applied

    # Set the threshold field as the background mesh field
    minimum = model.mesh.field.add("Min")
    model.mesh.field.setNumbers(minimum, "FieldsList", [threshold_field])
    model.mesh.field.setAsBackgroundMesh(threshold_field)
    # set defaults to zero to override them with the above settings
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

end

"""


"""
function _set_prog_suffix(proglc, lcmin, lcmax, distmin, distmax)
    if proglc
        suffix = "prog_lc_range_$(lcmin)_$(lcmax)_dist_range_$(distmin)_$(distmax)_"
    else
        suffix = ""
    end
    suffix
end
