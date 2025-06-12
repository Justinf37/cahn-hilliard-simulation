using Gmsh
using Gridap, GridapGmsh
using WriteVTK


Gmsh.initialize()
gmsh.model.add("box.geo")


vol = gmsh.model.occ.addBox(0.0,0.0,0.0,1.0,1.0,0.2)
gmsh.model.occ.synchronize()
model = gmsh.model
model.addPhysicalGroup(2, [-1], 1, "heated")
model.addPhysicalGroup(3, [vol],      2, "volume")
h=0.02

gmsh.model.mesh.generate(3)
gmsh.model.mesh.refine()
gmsh.model.mesh.refine()

gmsh.write("demo/box.msh")
Gmsh.finalize()


order = 1
model = GmshDiscreteModel("demo/box.msh")
reffe = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(
  model, reffe;
  dirichlet_tags = "heated"             
)


Th = 1.0
g(t) = x -> Th
U    = TransientTrialFESpace(V, g)


Ω  = Triangulation(model)
dΩ = Measure(Ω, 2*order)

D = 1.0 
M(t,dtu,v) = ∫( dtu * v )dΩ
K(t,u,v)   = ∫( D * ∇(u) ⋅ ∇(v) )dΩ
l(t,v)      = ∫( 0.0 * v )dΩ

op = TransientLinearFEOperator((K,M), l, U, V)

u0(x) = 0.0
uh0  = interpolate_everywhere(u0, U(0.0))


Δt = 0.01; θ = 0.5; t0 = 0.0; tF = 1.0
solver = ThetaMethod(LUSolver(), Δt, θ)
uh = solve(solver, op, t0, tF, uh0)


if !isdir("resultats/chaleur3D")
  mkdir("resultats/chaleur3D")
end


createpvd("resultats/chaleur3D/chaleur3D") do pvd

  pvd[0] = createvtk(Ω, "resultats/chaleur3D/chaleur3D_0.vtu",
                     cellfields = ["u" => uh0],
                     compress   = false,
                     append     = false)
 
  for (tn, uhn) in uh
    pvd[tn] = createvtk(Ω, "resultats/chaleur3D/chaleur3D_$tn.vtu",
                        cellfields = ["u" => uhn],
                        compress   = false,
                        append     = false)
  end
end

println("Simulation sur maillage Gmsh terminée ")
