using ForwardDiff
using WriteVTK
using Gridap
using Random

using Gridap.FESpaces: get_free_dof_values

Random.seed!(1234)


domain    = (0.0,1.0,0.0,1.0)
partition = (200,200)
model     = CartesianDiscreteModel(domain, partition)

order = 1
reffe  = ReferenceFE(lagrangian, Float64, order)

Vc = TestFESpace(model, reffe; conformity=:H1)  
Vμ = TestFESpace(model, reffe; conformity=:H1)  
Uc = TrialFESpace(Vc)                           
Uμ = TrialFESpace(Vμ)                           

X = MultiFieldFESpace([Uc, Uμ]) 
Y = MultiFieldFESpace([Vc, Vμ])  


Ω  = Triangulation(model)
dΩ = Measure(Ω, 2*order)

M = 1.0                      
λ     = 1e-2                    
f(c)   = 100*c^2*(1 - c)^2       
dfdc(c) = 100 * c * (1 - c) * (1 - 2c)


t0    = 0.0
tF    = 1e-3

c0fun(x) = 0.63 + 0.02*(rand() - 0.5)

uh₀_c = interpolate_everywhere(c0fun, Uc)
uh₀_μ = interpolate_everywhere(x -> 0.0, Uμ)
u0 = FEFunction(X, vcat(get_free_dof_values(uh₀_c), get_free_dof_values(uh₀_μ)))




m(t, (c,μ), (v,q)) = ∫( v*c )dΩ

a(t, (c,μ), (v,q)) = ∫(
  M * ∇(v) ⋅ ∇(μ) +
  q * μ -
  q * dfdc(c) -
  λ * ∇(q) ⋅ ∇(c)
)dΩ


ℓ(t, (v,q)) = ∫( 0.0*v + 0.0*q )dΩ

op = TransientLinearFEOperator(
  (a, m),  
  ℓ,      
  X, Y     
)
 
Δt = 5e-6
t0 = 0.0
tF = 1e-3
θ = 1.0
ls = LUSolver()
solver = ThetaMethod(ls, Δt, θ)
 
uh = solve(solver, op, t0, tF, u0)
ch,muh = uh

if !isdir("resultats/cahnhilliard")
  mkdir("resultats/cahnhilliard")
end


createpvd("resultats/cahnhilliard/res") do pvd


  pvd[0] = createvtk(Ω, "resultats/cahnhilliard/res_0.vtu",
                     cellfields = ["c" => u0[1], "mu" => u0[2]],
                     compress   = false,
                     append     = false)

 
  for (tn, uhn) in uh

    chn, muhn = uhn

    pvd[tn] = createvtk(Ω, "resultats/cahnhilliard/res_$tn.vtu",
                        cellfields = ["c" => chn, "mu" => muhn],
                        compress   = false,
                        append     = false)
  end

end
