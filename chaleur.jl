using Gmsh
using Gridap, GridapGmsh
using WriteVTK


Gmsh.initialize()
gmsh.model.add("rectangle.geo")


vol = gmsh.model.occ.addRectangle(0.0,0.0,0.0,1.0,1.0)
gmsh.model.occ.synchronize()
model = gmsh.model
model.addPhysicalGroup(1, [4], 1, "heated")
model.addPhysicalGroup(2, [vol],      1, "volume")
h=0.02

gmsh.model.mesh.generate(2)
gmsh.model.mesh.refine()
gmsh.model.mesh.refine()

gmsh.write("demo/rectangle.msh")
Gmsh.finalize()


order = 1
model = GmshDiscreteModel("demo/rectangle.msh")
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


D = 0.1   
M(t,dtu,v) = ∫( dtu * v )dΩ
K(t,u,v)   = ∫( D * ∇(u) ⋅ ∇(v) )dΩ
l(t,v)      = ∫( 0.0 * v )dΩ

op = TransientLinearFEOperator((K,M), l, U, V)


u0(x) = 0.0
uh0  = interpolate_everywhere(u0, U(0.0))


Δt = 0.01; θ = 0.5; t0 = 0.0; tF = 1.0
solver = ThetaMethod(LUSolver(), Δt, θ)
uh = solve(solver, op, t0, tF, uh0)


if !isdir("resultats/chaleur")
  mkdir("resultats/chaleur")
end


createpvd("resultats/chaleur/chaleur") do pvd
  
  pvd[0] = createvtk(Ω, "resultats/chaleur/chaleur_0.vtu",
                     cellfields = ["u" => uh0],
                     compress   = false,
                     append     = false)
  
  for (tn, uhn) in uh
    pvd[tn] = createvtk(Ω, "resultats/chaleur/chaleur_$tn.vtu",
                        cellfields = ["u" => uhn],
                        compress   = false,
                        append     = false)
  end
end

println("Simulation sur maillage Gmsh terminée ")
