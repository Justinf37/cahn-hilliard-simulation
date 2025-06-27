using ForwardDiff
using WriteVTK
using Gridap
using GridapODEs
using Random
import Sundials as snd
using Gridap.FESpaces: get_free_dof_values

Random.seed!(1234)


domain    = (0.0,1.0,0.0,1.0)
partition = (50,50)
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

D(x)   = TensorValue(1.0, 0.0, 0.0, 0.0)
D_field= CellField(D, Ω)


m(t, (c,μ), (v,q)) = ∫( v*c )dΩ

a(t, (c,μ), (v,q)) = ∫(
  ∇(v) ⋅ (D_field ⋅ ∇(μ))  +
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
function prep_solver_objects(op, init_fun, 
    tspan::Tuple{Float64, Float64}; 
    θ::Float64=0.1, dt::Float64=0.1)


    dt0 = θ * dt


    init_fun_snd = get_free_dof_values(init_fun)

  
    res!, jac! = GridapODEs.ODETools.diffeq_wrappers(op, init_fun_snd, tspan)

   
    J_ = GridapODEs.TransientFETools.prototype_jacobian(op, init_fun_snd, tspan)

   
    R_ = copy(init_fun_snd)

   
    res!(R_, init_fun_snd, init_fun_snd, [], tspan[1])
    jac!(J_, init_fun_snd, init_fun_snd, [], (1.0 / dt0), tspan[1])

   
    params = Dict(
        :θ     => θ,
        :dt    => dt,
        :dt0   => dt0,
        :tspan => tspan,
    )

    return R_, J_, res!, jac!, init_fun_snd, params
end


t0 = 0.0
tF = 1e-2
Δt  = 5e-6
θ   = 1.0        
nsteps = Int(round((tF - t0)/Δt)) + 1


R_, J_, res!, jac!, init_fun_snd, params =
  prep_solver_objects(op, u0, (t0, tF);
                      θ=θ, dt=Δt)


diff_vars = fill(true, length(init_fun_snd))


f_sundials = snd.DAEFunction{true,true}(res!;
                                        jac_prototype=J_,
                                        jac=jac!)
prob = snd.DAEProblem{true}(f_sundials,
                            init_fun_snd,     
                            init_fun_snd,      
                            (t0,tF),
                            (); differential_vars=diff_vars)


sol = snd.solve(prob,
                snd.IDA(linear_solver=:KLU, init_all=false);
                abstol=1e-12,
                reltol=1e-10,
                saveat=range(t0, tF, length=nsteps),
                progress=true)
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
