
using Gmsh
using Gridap, GridapGmsh
using WriteVTK
using Printf






Gmsh.initialize()
gmsh.model.add("rectangle.geo")


vol = gmsh.model.occ.addRectangle(0.0,0.0,0.0,1.0,1.0)
gmsh.model.occ.synchronize()
model = gmsh.model
model.addPhysicalGroup(1, [4], 1, "left")
model.addPhysicalGroup(1, [2], 2, "right")
model.addPhysicalGroup(2, [vol],      1, "volume")
h=0.02

gmsh.model.mesh.generate(2)


gmsh.write("demo/carre.msh")
Gmsh.finalize()


model = GmshDiscreteModel("demo/carre.msh")

order   = 1
reffe  = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(model, reffe; dirichlet_tags = ["left"])

U       = TransientTrialFESpace(V, t -> (x -> 0.0))


Ω  = Triangulation(model)
dΩ = Measure(Ω, 2*order)
Γ  = BoundaryTriangulation(model, tags=1)
dΓ = Measure(Γ, 2*order)

diffusivity = 0.001  
dt          = 0.05    
θ           = 0.5     
t0, tF      = 0.0, 10.0


k0    = 1e-2
α     = 0.7
nF_RT = 1.0
c_ref = 1.0
η     = 0.5

kBV = k0 * exp(α * nF_RT * η)
bBV = -k0 * c_ref * exp(-(1 - α) * nF_RT * η)


u0(x) = exp(-100 * ((x[1] - 0.5)^2 + (x[2] - 0.5)^2))
uh0   = interpolate_everywhere(u0, U(t0))


mass_form   = (t, du, v) -> ∫(du * v)dΩ
stiff_form  = (t, u, v)  -> ∫(diffusivity * ∇(u)⋅∇(v))dΩ + ∫(kBV * u * v)dΓ
source_form = (t, v)     -> ∫(bBV * v)dΓ


op = TransientLinearFEOperator((stiff_form, mass_form),
                               source_form,
                               U, V)


solver   = ThetaMethod(LUSolver(), dt, θ)
solution = solve(solver, op, t0, tF, uh0)


dir = "resultats/diffusionv3"
isdir(dir) || mkdir(dir)
createpvd(dir * "/res") do pvd
  pvd[0] = createvtk(Ω, dir * "/res_0000.vtu";
                     cellfields=["u"=>uh0], compress=false, append=false)
  for (time, uh) in solution
    file = @sprintf("%s/res_%06.2f.vtu", dir, time)
    pvd[time] = createvtk(Ω, file;
                          cellfields=["u"=>uh], compress=false, append=false)
  end
end
println("VTK export terminé dans ", dir)

