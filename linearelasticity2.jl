using Pkg
using Gmsh
using Gridap, GridapGmsh
using LinearAlgebra, WriteVTK

Gmsh.initialize()
gmsh.model.add("sphere_model")

R = 1.0
lc = 0.1

vol = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, R)
gmsh.model.occ.synchronize()

ctr = gmsh.model.occ.addPoint(0.0, 0.0, 0.0)
gmsh.model.occ.synchronize()

gmsh.model.mesh.embed(0, [ctr], 3, vol)

model = gmsh.model
model.addPhysicalGroup(0, [ctr],      1, "centre")
model.addPhysicalGroup(2, [-1],       1, "boundary1")
model.addPhysicalGroup(3, [vol],      1, "volume")

gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc)
gmsh.model.mesh.generate(3)
gmsh.write("demo/sphere_boundary_center.msh")
Gmsh.finalize()

model = GmshDiscreteModel("sphere_boundary_center.msh") 
order = 1
reffe = ReferenceFE(lagrangian, VectorValue{3,Float64}, order)

V = TestFESpace(model, reffe;
  conformity = :H1,
  dirichlet_tags= ["centre"],
  dirichlet_masks = [(true,true,true)]
)


U = TrialFESpace(V)


degree = 2 * order
Ω = Triangulation(model); dΩ = Measure(Ω, degree)
Γ = BoundaryTriangulation(model, tags = ["boundary1"])
dΓ = Measure(Γ, degree)


E = 70e9
ν = 0.33
λ = (E * ν) / ((1 + ν)*(1 - 2ν))
μ = E / (2 * (1 + ν))
I3 = one(TensorValue{3,3,Float64})
ε(u) = 1/2 * (∇(u) + ∇(u)')
σ(εu) = λ * tr(εu) * I3 + 2μ * εu


eps_stab = 0.1
a(u,v) = ∫( inner(ε(v), σ(ε(u))) + eps_stab * inner(∇(u), ∇(v)) )dΩ


t0 = 1e4
nΓ = get_normal_vector(Γ)
traction_exact(x) = t0 * (x / norm(x))
l(v) = ∫( t0 * inner(nΓ, v) )dΓ


op = AffineFEOperator(a, l, U, V)
uh = solve(op)




u_analytic(x) = begin
  r = norm(x)
  if r < 1e-12
    return VectorValue(0.0, 0.0, 0.0)
  end
  factor = (1 + ν)*(1 - 2ν) / (E * (1 - ν))
  u_r = t0 * r * factor
  return u_r * (x / r)
end


uh_exact = interpolate_everywhere(u_analytic, V)



e_field = uh - uh_exact







writevtk(Ω, "results/compare_uh";
  cellfields = [
    "uh_numeric" => uh,
    "uh_exact" => uh_exact,
    "error" => e_field
  ],
  compress = false,
        append   = false)

