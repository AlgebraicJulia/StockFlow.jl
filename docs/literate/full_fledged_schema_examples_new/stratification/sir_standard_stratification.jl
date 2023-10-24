using GraphViz

using StockFlow
using StockFlow.Syntax

using Catlab
using Catlab.CategoricalAlgebra
using LabelledArrays
using OrdinaryDiffEq
using Plots

using Catlab.Graphics
using Catlab.Programs
using Catlab.Theories
using Catlab.WiringDiagrams

using Catlab.Graphics.Graphviz: Html
using Catlab.Graphics.Graphviz



# Functions for graphing typed Petri nets
colors_vflow = ["antiquewhite4","antiquewhite", "gold", "saddlebrown", "slateblue", "blueviolet", "olive"]
colors_s = ["deeppink","darkorchid","darkred","coral"] # red series
colors_sv = ["cornflowerblue","cyan4","cyan","chartreuse"] # green and blue series
colors_p = ["gold","gold4","darkorange1","lightgoldenrod","goldenrod"] # yellow and orange

flatten(fname::Symbol) = "$fname"

function flatten(fname::Tuple)
    names = split(replace(string(fname), "("=>"", ")"=>"", ":"=>""), ",")
    for i in 1:length(names)
        name = strip(names[i])
        if name[1:2] == "id"
            continue
        end
        return name
    end
    return "id"
end

def_stock(typed_StockFlow::ACSetTransformation, colors) = 
  (p,s) -> ("s$s", Attributes(:label=>sname(p,s) isa Tuple where T ? Html(replace(string(sname(p,s)), ":"=>"", "," => "<BR/>", "("=>"", ")"=>"")) : "$(sname(p,s))",
                                     :shape=>"square", 
                                     :color=>"black", 
                                     :style=>"filled", 
                                     :fillcolor=>colors[typed_StockFlow[:S](s)]))

def_parameter(typed_StockFlow::ACSetTransformation, colors) = 
(p, pp) -> ("p$pp", Attributes(:label=>pname(p,pp) isa Tuple where T ? Html(replace(string(pname(p,pp)), ":"=>"", "," => "<BR/>", "("=>"", ")"=>"")) : "$(pname(p,pp))",
                                     :shape=>"circle", 
                                     :color=>colors[typed_StockFlow[:P](pp)],
                                     :fontcolor=>colors[typed_StockFlow[:P](pp)]))

def_auxiliaryVF(typed_StockFlow::ACSetTransformation, colors)=
  (p, v) -> ("v$v", Attributes(:label=>make_v_expr(p,v) isa Tuple where T ? Html(replace(string(make_v_expr(p,v)), ":"=>"", "," => "<BR/>", "("=>"", ")"=>"")) : "$(make_v_expr(p,v))",
                                          :shape=>"plaintext", 
                                          :fontcolor=>colors[typed_StockFlow[:V](v)]))


def_sumV(typed_StockFlow::ACSetTransformation, colors) = 
  (p, sv) -> ("sv$sv", Attributes(:label=>svname(p,sv) isa Tuple where T ? Html(replace(string(svname(p,sv)), ":"=>"", "," => "<BR/>", "("=>"", ")"=>"")) : "$(svname(p,sv))",
                                       :shape=>"circle", 
                                       :color=>"black",
                                       :fillcolor=>colors[typed_StockFlow[:SV](sv)], 
                                       :style=>"filled"))  


def_flow_V(typed_StockFlow::ACSetTransformation, colors)=
  (p, us, ds, v, f) -> begin
    labelfontsize = "6"
    colorType = colors[typed_StockFlow[:F](f)]
    color = "$colorType"*":invis:"*"$colorType"
    arrowhead = "none"
    splines = "ortho"
    return ([us, "v$v"],Attributes(:label=>"", :labelfontsize=>labelfontsize, :color=>color, :arrowhead=>arrowhead, :splines=>splines)),
           (["v$v", ds],Attributes(:label=>Html(flatten(fname(p,f))), :labelfontsize=>labelfontsize, :color=>color, :splines=>splines))            
end
        
def_flow_noneV(typed_StockFlow::ACSetTransformation, colors)=
  (p, us, ds, f) -> begin
     colorType = colors[typed_StockFlow[:F](f)]
     color = "$colorType"*":invis:"*"$colorType"   
     ([us, ds],Attributes(:label=>Html(flatten(fname(p,f))), :labelfontsize=>"6", :color=>color))
end

def_flow_V(colors = colors_vflow)=
  (p, us, ds, v, f) -> begin
    labelfontsize = "6"
    colorType = colors[f]
    color = "$colorType"*":invis:"*"$colorType"
    arrowhead = "none"
    splines = "ortho"
    return ([us, "v$v"],Attributes(:label=>"", :labelfontsize=>labelfontsize, :color=>color, :arrowhead=>arrowhead, :splines=>splines)),
           (["v$v", ds],Attributes(:label=>Html(flatten(fname(p,f))), :labelfontsize=>labelfontsize, :color=>color, :splines=>splines))            
    
end
        
def_flow_noneV(colors = colors_vflow)=
  (p, us, ds, f) -> begin
     colorType = colors[f]
     color = "$colorType"*":invis:"*"$colorType"   
     ([us, ds],Attributes(:label=>Html(flatten(fname(p,f))), :labelfontsize=>"6", :color=>color))
end

GraphF_typed(typed_StockFlow::ACSetTransformation, colors_vflow = colors_vflow, colors_s = colors_s, colors_p = colors_p, colors_sv = colors_sv; schema::String="C", type::String="SFVL", rd::String="LR") = GraphF(dom(typed_StockFlow),
    make_stock = def_stock(typed_StockFlow, colors_s), make_auxiliaryV=def_auxiliaryVF(typed_StockFlow, colors_vflow), make_sumV=def_sumV(typed_StockFlow, colors_sv), 
    make_flow_V=def_flow_V(typed_StockFlow, colors_vflow), make_flow_noneV=def_flow_noneV(typed_StockFlow, colors_vflow),make_parameter=def_parameter(typed_StockFlow, colors_p),schema=schema, type=type, rd=rd
)

# using StockFlow: def_stock, def_parameter, def_auxiliaryVF, def_sumV, def_flow_V

# # Functions for graphing typed Petri nets
# colors_vflow = ["antiquewhite4","antiquewhite", "gold", "saddlebrown", "slateblue", "blueviolet", "olive"]
# colors_s = ["deeppink","darkorchid","darkred","coral"] # red series
# colors_sv = ["cornflowerblue","cyan4","cyan","chartreuse"] # green and blue series
# colors_p = ["gold","gold4","darkorange1","lightgoldenrod","goldenrod"] # yellow and orange


# flatten(fname::Symbol) = "$fname"

# function flatten(fname::Tuple)
#     names = split(replace(string(fname), "("=>"", ")"=>"", ":"=>""), ",")
#     for i in 1:length(names)
#         name = strip(names[i])
#         if name[1:2] == "id"
#             continue
#         end
#         return name
#     end
#     return "id"
# end


# def_flow_noneV(colors = colors_vflow)=
#   (p, us, ds, f) -> begin
#      colorType = colors[f]
#      color = "$colorType"*":invis:"*"$colorType"   
#      ([us, ds],Attributes(:label=>Html(flatten(fname(p,f))), :labelfontsize=>"6", :color=>color))
# end

# GraphF_typed(typed_StockFlow::ACSetTransformation, colors_vflow = colors_vflow, colors_s = colors_s, colors_p = colors_p, colors_sv = colors_sv; schema::String="C", type::String="SFVL", rd::String="LR") = GraphF(dom(typed_StockFlow),
#     make_stock = def_stock(typed_StockFlow, colors_s), make_auxiliaryV=def_auxiliaryVF(typed_StockFlow, colors_vflow), make_sumV=def_sumV(typed_StockFlow, colors_sv), 
#     make_flow_V=def_flow_V(typed_StockFlow, colors_vflow), make_flow_noneV=def_flow_noneV(typed_StockFlow, colors_vflow),make_parameter=def_parameter(typed_StockFlow, colors_p),schema=schema, type=type, rd=rd
# )


s_type = @stock_and_flow begin
    :stocks
    pop
    
    :parameters
    c
    β
    rFstOrder
    rAge
    
    :dynamic_variables
    v_prevalence = pop / N
    v_meanInfectiousContactsPerS = c * v_prevalence
    v_perSIncidenceRate = β * v_meanInfectiousContactsPerS
    v_inf = pop * v_perSIncidenceRate
    v_fstOrder = pop * rFstOrder
    v_aging = pop * rAge
    
    :flows
    pop => f_inf(v_inf) => pop
    pop => f_fstOrder(v_fstOrder) => pop
    pop => f_aging(v_aging) => pop

    
    :sums
    N = [pop]
    
    
end


GraphF(s_type)

GraphF_typed(id(s_type))

# eliminate the attribute of name to enable pass the natural check
# only eliminate the name, the other two attributes should be okay
s_type = map(s_type, Name=name->nothing, Op=op->nothing, Position=pos->nothing);

s, = parts(s_type, :S)
N, = parts(s_type, :SV)
lsn, = parts(s_type, :LS)
f_inf, f_fstorder, f_aging = parts(s_type, :F)
i_inf, i_fstorder, i_aging = parts(s_type, :I) # note, different order from previous in both inflow and outflow
o_inf, o_fstorder, o_aging = parts(s_type, :O)
v_IN, v_cIN, v_betacIN, v_inf, v_fstOrder, v_aging = parts(s_type, :V)
lv_IN1, lv_inf1, lv_fstOrder1, lv_aging1 = parts(s_type, :LV)
lsv_IN2, = parts(s_type, :LSV)
p_c, p_beta, p_rfstOrder, p_rAge = parts(s_type, :P)
lvv_cIN2, lvv_betacIN2, lvv_inf2 = parts(s_type, :LVV)
lpv_cIN1, lpv_betacIN1, lpv_fstOrder2, lpv_aging2 = parts(s_type, :LPV)

sir = @stock_and_flow begin
    :stocks
    S
    I
    R
    
    :parameters
    c
    β
    rRec
    rAge
    
    :dynamic_variables
    v_prevalence = I / N
    v_meanInfectiousContactsPerS = c * v_prevalence
    v_perSIncidenceRate = β * v_meanInfectiousContactsPerS
    v_newInfections = S * v_perSIncidenceRate
    v_newRecovery = I * rRec
    v_idS = S * rAge
    v_idI = I * rAge
    v_idR = R * rAge
    
    :flows
    S => f_idS(v_idS) => S
    S => f_inf(v_newInfections) => I
    I => f_idI(v_idI) => I
    I => f_rec(v_newRecovery) => R
    R => f_idR(v_idR) => R
    
    :sums
    N = [S, I, R]
    
    
end

GraphF(sir)

typed_aggregate_model=ACSetTransformation(sir, s_type,
  S = [s,s,s],
  SV = [N],
  LS = [lsn,lsn,lsn],   
  F = [f_aging, f_inf, f_aging, f_fstorder, f_aging],    
  I = [i_aging, i_inf, i_aging, i_fstorder, i_aging], #i_inf, i_fstorder, i_aging
  O = [o_aging, o_inf, o_aging, o_fstorder, o_aging],
  V = [v_IN, v_cIN, v_betacIN, v_inf, v_fstOrder, v_aging, v_aging, v_aging],
  LV = [lv_IN1, lv_inf1, lv_fstOrder1, lv_aging1, lv_aging1, lv_aging1],
  LSV = [lsv_IN2],
  P = [p_c, p_beta, p_rfstOrder, p_rAge],
  LVV = [lvv_cIN2, lvv_betacIN2, lvv_inf2],
  LPV = [lpv_cIN1, lpv_betacIN1, lpv_fstOrder2, lpv_aging2, lpv_aging2, lpv_aging2],
  Name = name -> nothing, Op=op->nothing, Position=pos->nothing
);
@assert is_natural(typed_aggregate_model)

GraphF_typed(typed_aggregate_model)

age2 = @stock_and_flow begin
    :stocks
    Child
    Adult
    
    :parameters
    c_C
    β
    r
    rAge
    c_A
    
    :dynamic_variables
    v_INC = Child / NC
    v_cINC = c_C * v_INC
    v_cβINC = β * v_cINC
    
    v_infC = Child * v_cβINC
    v_fstC = Child * r
    v_agingC = Child * rAge
    
    
    v_INA = Adult / NA
    v_cINA = c_A * v_INA
    v_cβINA = β * v_cINA
    
    v_infA = Adult * v_cβINA
    v_fstA = Adult * r
    
    :flows
    Child => f_infC(v_infC) => Child
    Child => f_frsC(v_fstC) => Child
    Child => f_aging(v_agingC) => Adult
    Adult => f_infA(v_infA) => Adult
    Adult => f_frsA(v_fstA) => Adult
    
    
    :sums
    NC = [Child]
    NA = [Adult]
end

typed_age_model=ACSetTransformation(age2, s_type,
  S = [s,s],
  SV = [N,N],
  LS = [lsn,lsn],   
  F = [f_inf, f_fstorder, f_aging, f_inf, f_fstorder],    
  I = [i_inf, i_fstorder, i_aging, i_inf, i_fstorder], 
  O = [o_inf, o_fstorder, o_aging, o_inf, o_fstorder],
  V = [v_IN, v_cIN, v_betacIN, v_inf, v_fstOrder, v_aging, v_IN, v_cIN, v_betacIN, v_inf, v_fstOrder],
  LV = [lv_IN1, lv_inf1, lv_fstOrder1, lv_aging1, lv_IN1, lv_inf1, lv_fstOrder1],
  LSV = [lsv_IN2, lsv_IN2],
  P = [p_c, p_beta, p_rfstOrder, p_rAge, p_c],
  LVV = [lvv_cIN2, lvv_betacIN2, lvv_inf2, lvv_cIN2, lvv_betacIN2, lvv_inf2],
  LPV = [lpv_cIN1, lpv_betacIN1, lpv_fstOrder2, lpv_aging2, lpv_cIN1, lpv_betacIN1, lpv_fstOrder2],
  Name = name -> nothing, Op=op->nothing, Position=pos->nothing
);
@assert is_natural(typed_age_model)

GraphF_typed(typed_age_model)

aged_sir = pullback(typed_aggregate_model, typed_age_model) |> apex |> rebuildStratifiedModelByFlattenSymbols;

GraphF(aged_sir)

LS = @stock_and_flow begin
    :stocks
    SChild
    IChild
    SAdult
    IAdult
    
    :parameters
    cc_C
    cc_A
    
    :dynamic_variables
    v_prevalencev_INC = IChild / NNC
    v_prevalencev_INA = IAdult / NNA
    v_meanInfectiousContactsPerSv_cINC = cc_C * v_prevalencev_INC
    v_meanInfectiousContactsPerSv_cINA = cc_A * v_prevalencev_INA
    
    
    :sums
    NNC = [SChild, IChild]
    NNA = [SAdult, IAdult]
    
    
end

GraphF(LS)

IS = @stock_and_flow begin
    :stocks
    SChild
    IChild
    SAdult
    IAdult
    
    :parameters
    cc_C
    cc_A
    
    :dynamic_variables
    v_prevalencev_INC = IChild / NNC
    v_prevalencev_INA = IAdult / NNA
    v_meanInfectiousContactsPerSv_cINC = *(cc_C)
    v_meanInfectiousContactsPerSv_cINA = *(cc_A)
    
    :sums
    NNC = [SChild, IChild]
    NNA = [SAdult, IAdult]
    
    
end

GraphF(IS)

RS = @stock_and_flow begin
    :stocks
    SChild
    IChild
    SAdult
    IAdult
    
    :parameters
    fcc
    fca
    fac
    faa
    cc_C
    cc_A
    
    :dynamic_variables
    v_prevalencev_INC = IChild / NNC
    v_prevalencev_INA = IAdult / NNA
    v_CCContacts = fcc * v_prevalencev_INC
    v_CAContacts = fca * v_prevalencev_INA
    
    v_ACContacts = fac * v_prevalencev_INC
    v_AAContacts = faa * v_prevalencev_INA
    
    v_prevalencev_INC_post = v_CCContacts + v_CAContacts
    v_prevalencev_INA_post = v_ACContacts + v_AAContacts
    v_meanInfectiousContactsPerSv_cINC = cc_C * v_prevalencev_INC_post
    v_meanInfectiousContactsPerSv_cINA = cc_A * v_prevalencev_INA_post
    
    :sums
    NNC = [SChild, IChild]
    NNA = [SAdult, IAdult]
    
    
end


GraphF(RS)

using AlgebraicRewriting
using AlgebraicRewriting: rewrite
const hom = Catlab.CategoricalAlgebra.homomorphism

rule_S = Rule(hom(IS,LS), hom(IS,RS))
aged_sir_rewritten = rewrite(rule_S, aged_sir)
GraphF(aged_sir_rewritten)

# define values of constant parameters
p_stratified_sir = LVector(
    fcc=0.8, fca=0.2, fac=0.2, faa=0.8, cc_C=0.45, cc_A=0.55,
    ββ=0.8, rRecr=1.0/14.0, rAgerAge=1.0/(15.0*365.0) #which means the child age group is from 0 to 15 years old
)
# define initial values for stocks
u0_stratified_sir = LVector(
    SChild=990.0, IChild=10.0, RChild=0.0,
    SAdult=4900.0, IAdult=10.0, RAdult=0.0
);

prob_stratified_sir = ODEProblem(vectorfield(aged_sir_rewritten),u0_stratified_sir,(0.0,100.0),p_stratified_sir);
sol_stratified_sir = solve(prob_stratified_sir,Tsit5(),abstol=1e-8);
plot(sol_stratified_sir)

# to have the figures plotted fix to the wider of the cells
HTML("""
<style>
.output_svg div{
  width: 100% !important;
  height: 100% !important;
}
</style>
""")



