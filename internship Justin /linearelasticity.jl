
using Gridap, GridapGmsh, WriteVTK




model = GmshDiscreteModel("demo/cercle4.msh")


order = 1
reffe = ReferenceFE(lagrangian, VectorValue{3,Float64}, order)
V     = TestFESpace(
            model, reffe;
            conformity      = :H1,
            dirichlet_tags  = ["boundary1", "boundary2"],
            dirichlet_masks = [(true,false,false), (true,true,true)]
        )


g1(x) = VectorValue(0.5, 0.0, 0.0) 
g2(x) = VectorValue(0.0,   0.0, 0.0)  
U = TrialFESpace(V, [g1, g2])


const E  = 70.0e9      
const ν  = 0.33       
const λ  = (E*ν)/((1+ν)*(1-2ν))
const μ  = E/(2*(1+ν))




const I3 = one(TensorValue{3,3,Float64})

ε(u) = 1/2*(gradient(u) + transpose(gradient(u)))
σ(ε) = λ*tr(ε)*I3 + 2μ*ε


degree = 2*order
Ω      = Triangulation(model)
dΩ     = Measure(Ω, degree)

a(u,v) = ∫( ε(v) ⊙ σ(ε(u)) )dΩ
l(v)    = 0 


op = AffineFEOperator(a, l, U, V)
uh = solve(op)


writevtk(Ω, "resultats/cercleimpact";
        cellfields = [
          "uh"    => uh,                    
          "ε(uh)" => ε(uh),                  
          "σ(uh)" => σ(ε(uh))                
        ],
        compress = false,
        append   = false)
