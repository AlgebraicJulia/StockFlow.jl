# # SIR Linear Stratification

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

l_type = @stock_and_flow begin 
    :stocks
    pop
    
    :parameters
    μ
    δ
    rFstOrder
    rage
    
    :dynamic_variables
    v_aging = pop * rage
    v_fstOrder = pop * rFstOrder
    v_birth = N * μ
    v_death = pop * δ
    
    :flows
    pop => f_aging(v_aging) => pop
    pop => f_fstOrder(v_fstOrder) => pop
    CLOUD => f_birth(v_birth) => pop
    pop => f_death(v_death) => CLOUD
    
    :sums
    N = [pop]
    
end

GraphF_typed(id(l_type))

# eliminate the attribute of name to enable pass the natural check
# only eliminate the name, the other two attributes should be okay
l_type = map(l_type, Name=name->nothing, Op=op->nothing, Position=pos->nothing);

s, = parts(l_type, :S)
N, = parts(l_type, :SV)
lsn, = parts(l_type, :LS)
f_aging, f_fstorder, f_birth, f_death = parts(l_type, :F)
i_aging, i_fstorder, i_birth = parts(l_type, :I)
o_aging, o_fstorder, o_death = parts(l_type, :O)
v_aging, v_fstorder, v_birth, v_death = parts(l_type, :V)
lv_aging1, lv_fstorder1, lv_death1 = parts(l_type, :LV)
lsv_birth1, = parts(l_type, :LSV)
p_μ, p_δ, p_rfstOrder, p_rage = parts(l_type, :P)
lpv_aging2, lpv_fstorder2, lpv_birth2, lpv_death2 = parts(l_type, :LPV)


WeightModel = @stock_and_flow begin
    :stocks
    NormalWeight
    OverWeight
    Obese
    
    :parameters
    μ
    δw
    rw
    ro
    δo
    rage
    
    :dynamic_variables
    v_NewBorn = N * μ
    v_DeathNormalWeight = NormalWeight * δw
    v_BecomingOverWeight = NormalWeight * rw
    v_DeathOverWeight = OverWeight * δw
    v_BecomingObese = OverWeight * ro
    v_DeathObese = Obese * δo
    v_idNW = NormalWeight * rage
    v_idOW = OverWeight * rage
    v_idOb = Obese * rage
    
    :flows
    CLOUD => f_NewBorn(v_NewBorn) => NormalWeight
    NormalWeight => f_DeathNormalWeight(v_DeathNormalWeight) => CLOUD
    NormalWeight => f_BecomingOverWeight(v_BecomingOverWeight) => OverWeight
    OverWeight => f_DeathOverWeight(v_DeathOverWeight) => CLOUD
    
    OverWeight => f_BecomingObese(v_BecomingObese) => Obese
    Obese => f_DeathObese(v_DeathObese) => CLOUD
    NormalWeight => f_idNW(v_idNW) => NormalWeight
    OverWeight => f_idOW(v_idOW) => OverWeight
    Obese => f_idOb(v_idOb) => Obese
    
    :sums
    N = [NormalWeight, OverWeight, Obese]
    
end

GraphF(WeightModel, rd="TB")

typed_WeightModel=ACSetTransformation(WeightModel, l_type,
  S = [s,s,s],
  SV = [N],
  LS = [lsn,lsn,lsn],   
  F = [f_birth, f_death, f_fstorder, f_death, f_fstorder, f_death, f_aging, f_aging, f_aging],    
  I = [i_birth, i_aging, i_fstorder, i_aging, i_fstorder, i_aging], 
  O = [o_death, o_fstorder, o_aging, o_death, o_fstorder, o_aging, o_death, o_aging],
  V = [v_birth, v_death, v_fstorder, v_death, v_fstorder, v_death, v_aging, v_aging, v_aging],
  LV = [lv_death1, lv_fstorder1, lv_death1, lv_fstorder1, lv_death1, lv_aging1, lv_aging1, lv_aging1],
  LSV = [lsv_birth1],
  P = [p_μ, p_δ, p_rfstOrder, p_rfstOrder, p_δ, p_rage],
  LPV = [lpv_birth2, lpv_death2, lpv_fstorder2, lpv_death2, lpv_fstorder2, lpv_death2, lpv_aging2, lpv_aging2, lpv_aging2],
  Name = name -> nothing, Op=op->nothing, Position=pos->nothing
);
@assert is_natural(typed_WeightModel)
GraphF_typed(typed_WeightModel, rd="TB")

ageWeightModel = @stock_and_flow begin
    :stocks
    Child
    Adult
    Senior
    
    :parameters
    μ
    δC
    δA
    δS
    rageCA
    rageAS
    r
    
    :dynamic_variables
    v_NB = N * μ
    v_DeathC = Child * δC
    v_idC = Child * r
    v_agingCA = Child * rageCA
    v_DeathA = Adult * δA
    v_idA = Adult * r
    v_agingAS = Adult * rageAS
    v_DeathS = Senior * δS
    v_idS = Senior * r
    
    :flows
    CLOUD => f_NB(v_NB) => Child
    Child => f_idC(v_idC) => Child
    Child => f_DeathC(v_DeathC) => CLOUD
    Child => f_agingCA(v_agingCA) => Adult
    Adult => f_idA(v_idA) => Adult
    Adult => f_DeathA(v_DeathA) => CLOUD
    Adult => f_agingAS(v_agingAS) => Senior
    Senior => f_idS(v_idS) => Senior
    Senior => f_DeathS(v_DeathS) => CLOUD
    
    :sums
    N = [Child, Adult, Senior]
    
    
end


GraphF(ageWeightModel)

typed_ageWeightModel=ACSetTransformation(ageWeightModel, l_type,
  S = [s,s,s],
  SV = [N],
  LS = [lsn,lsn,lsn],   
  F = [f_birth, f_fstorder, f_death, f_aging, f_fstorder, f_death, f_aging, f_fstorder, f_death],    
  I = [i_birth, i_fstorder, i_aging, i_fstorder, i_aging, i_fstorder], 
O = [o_fstorder, o_death, o_aging, o_fstorder, o_death, o_aging, o_fstorder, o_death],
V = [v_birth, v_death, v_fstorder, v_aging, v_death, v_fstorder, v_aging, v_death, v_fstorder],
  LV = [lv_death1, lv_fstorder1, lv_aging1, lv_death1, lv_fstorder1, lv_aging1, lv_death1, lv_fstorder1],
  LSV = [lsv_birth1],
  P = [p_μ, p_δ, p_δ, p_δ, p_rage, p_rage, p_rfstOrder],
  LPV = [lpv_birth2, lpv_death2, lpv_fstorder2, lpv_aging2, lpv_death2, lpv_fstorder2, lpv_aging2, lpv_death2, lpv_fstorder2],
  Name = name -> nothing, Op=op->nothing, Position=pos->nothing
);
@assert is_natural(typed_ageWeightModel)
GraphF_typed(typed_ageWeightModel)

aged_weight = pullback(typed_WeightModel, typed_ageWeightModel) |> apex |> rebuildStratifiedModelByFlattenSymbols;
GraphF(aged_weight)

p_weight = LVector(
    μμ=12.5/1000,δwδC=2.0/1000,δoδC=8.0/1000,δwδA=4.0/1000,δoδA=13.0/1000,δwδS=8.0/1000,δoδS=30.0/1000,
    ragerageCA=1.0/(12.0*365.0),ragerageAS=1.0/(50.0*365.0),rwr=0.03,ror=0.06
)

u0_weight = LVector(
    NormalWeightChild=95811.0*12.0/82.0, OverWeightChild=27709.0*12.0/82.0, ObeseChild=30770.0*12.0/82.0,
    NormalWeightAdult=95811.0*50.0/82.0, OverWeightAdult=27709.0*50.0/82.0, ObeseAdult=30770.0*50.0/82.0,
    NormalWeightSenior=95811.0*20.0/82.0, OverWeightSenior=27709.0*20.0/82.0, ObeseSenior=30770.0*20.0/82.0
);

prob_stratified_weight = ODEProblem(vectorfield(aged_weight),u0_weight,(0.0,100.0),p_weight);
sol_stratified_weight = solve(prob_stratified_weight,Tsit5(),abstol=1e-8);
plot(sol_stratified_weight)

## to have the figures plotted fix to the wider of the cells
## HTML("""
## <style>
## .output_svg div{
##   width: 100% !important;
##   height: 100% !important;
## }
## </style>
## """)



