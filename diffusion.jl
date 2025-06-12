using Gridap
using WriteVTK


domain    = (0.0, 1.0, 0.0, 1.0)
partition = (20, 20)
model     = CartesianDiscreteModel(domain, partition)


order = 1
reffe = ReferenceFE(lagrangian, Float64, order)


V = TestFESpace(model, reffe; dirichlet_tags="boundary")


g(t) = x -> 0.0  
U    = TransientTrialFESpace(V, g)


Ω  = Triangulation(model)
dΩ = Measure(Ω, 2*order)
Γ  = BoundaryTriangulation(model, tags="boundary")
dΓ = Measure(Γ, 2*order)


D   = 0.001        
Δt  = 0.05      
θ   = 0.5        
t0  = 0.0        
tF  = 10.0       


u0(x) = exp(-100*((x[1] - 0.5)^2 + (x[2] - 0.5)^2))
uh0 = interpolate_everywhere(u0, U(t0))


M(t, dtu, v) = ∫( dtu * v )dΩ
K(t, u, v)   = ∫( D * ∇(u) ⋅ ∇(v) )dΩ
l(t, v)      = ∫( 0.0 * v )dΩ  

op = TransientLinearFEOperator((K, M), l, U, V)



ls     = LUSolver()
solver = ThetaMethod(ls, Δt, θ)


uh = solve(solver, op, t0, tF, uh0)


if !isdir("resultats/diffusionv3")
  mkdir("resultats/diffusionv3")
end


createpvd("resultats/diffusionv3/res") do pvd

  pvd[0] = createvtk(Ω, "resultats/diffusionv3/res_0.vtu",
                     cellfields = ["u" => uh0],
                     compress   = false,
                     append     = false)

  for (tn, uhn) in uh
    pvd[tn] = createvtk(Ω, "resultats/diffusionv3/res_$tn.vtu",
                        cellfields = ["u" => uhn],
                        compress   = false,
                        append     = false)
  end
end

println("VTK export terminé dans le dossier diffusionv3/")


Nt = Int(round((tF - t0)/Δt))


u_ex(x, t) = (1 / (1 + 400*D*t)) * exp(
  -100 / (1 + 400*D*t) * ((x[1] - 0.5)^2 + (x[2] - 0.5)^2)
)


if !isdir("resultats/exact_solution")
  mkdir("resultats/exact_solution")
end

createpvd("resultats/exact_solution/exact") do pvd

  uex0 = interpolate_everywhere(x -> u_ex(x, t0), V)
  pvd[t0] = createvtk(Ω, "resultats/exact_solution/exact_0.vtu";
                      cellfields = ["u_exact" => uex0],
                      compress   = false,
                      append     = false)

  for i in 1:Nt
    tn   = t0 + i*Δt
    uexn = interpolate_everywhere(x -> u_ex(x, tn), V)
    pvd[tn] = createvtk(Ω, "resultats/exact_solution/exact_$(lpad(i,4,'0')).vtu";
                        cellfields = ["u_exact" => uexn],
                        compress   = false,
                        append     = false)
  end
end

println("Export de la solution exacte terminée")