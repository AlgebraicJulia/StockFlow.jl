# Domain Specific Languages

```julia
using StockFlow
using StockFlow.Syntax # stock_and_flow, foot, feet
using StockFlow.Syntax.Stratification # stratify, n_stratify
using StockFlow.Syntax.Composition # compose
using StockFlow.Syntax.Rewrite # rewrite
using StockFlow.Syntax.Homomorphism # hom
```

# Syntax

## @stock\_and\_flow

Creates a StockAndFlowF based upon the provided block expression.
Block has five headers you can use in it:

```julia
@stock_and_flow begin
        # to specify stock names
        :stocks
        S
        I
        R

        # to specify parameter names, which will be given values from outside
        # the model when evaluated
        :parameters
        β
        λ
        δ

        # to specify sum variables, which have at any given time have value
        # equal to the sum of all their linked stocks
        :sums
        N = [S, I, R]
        NI = [I]
        None = [] # 0

        # to specify how stocks, parameters, sums and other dynamic variables
        # relate to each other.  We use these to determine flow rates.

        # If we specify a dynamic variable with multiple operators, it'll
        # automatically be split into multiple.
        :dynamic_variables
        v₁ = S + I / N # this will be split into two dyvars
        v₂ = *(β, I)
        v₃ = v₂ + 1
        v₄ = log(S)

        # to specify the rate at which members of one population will
        # leave to become members of another population.
        # Use CLOUD or ☁ (\:cloud:) to indicate entering or leaving model.

        # You can also use the dynamic variable format inside a flow to create
        # a dynamic variable
        :flows
        CLOUD => f1(v1) => S
        S => Deaths_S(δ * S) => CLOUD
        S => infections(+λ) => I
        I => recoveries(v₄ + v₃ / v₂) => R

    end
```

The specific values in the above example are chosen to show the capabilities
of the syntax.  For more realistic examples, see src/PremadeModels.jl.

At present, any julia function is allowed to act as an operator for a 
dynamic variable, but in the future it'll be restricted to the below list:

Binary:  :+, :-, :*, :/, :÷, :^, :%, :log

Unary: :log, :exp, :sqrt, :+, :-
:plus\_one, :minus\_one, :reciprocal, :one\_minus, :plus\_two,
:minus_two

The unary functions can be expressed as plus\_one(X), or how you'd expect,
as X + 1 or 1 + X.

---

## @foot

Creates a StockAndFlow0 with stocks, sums, and links between them.
Takes the form @foot mystock => Y, B => Y, ..., where in each pair, the
first symbol is a stock and the second is a sum.  Multiple occurences
of the same symbol will be interpreted as the same instance of that sum
or stock.  If you have the same link multiple times, it will not be
deduplicated.  Use () => Y to indicate no stock, X => () to indicate no sum,
or () => () to indicate an empty foot.

```julia
(@foot () => ()) == StockAndFlow0([],[],[])
(@foot () => X) == StockAndFlow0([],[:X],[])
(@foot Y => ()) == StockAndFlow0([:Y],[],[])
(@foot Y => X) == StockAndFlow0([:Y],[:X],[:Y => :X])

(@foot Y => X, () => ()) == StockAndFlow0([:Y],[:X],[:Y => :X])
(@foot Y => X, Y => X) == StockAndFlow0([:Y],[:X],[:Y => :X, :Y => :X])
(@foot Y => X, A => X) == StockAndFlow0([:Y, :A],[:X],[:Y => :X, :A => :X])
(@foot Y => X, Y => Z) == StockAndFlow0([:Y],[:X, :Z],[:Y => :X, :Y => :Z])
```

---

## @feet

Given a block expression, produces an array of feet for each line.
An empty block produces an empty array.

```julia
(@feet begin end) == Vector{StockAndFlow0}()

(@feet begin 
() => ()
end) == [StockAndFlow0([],[],[])]

(@feet begin 
A => B, AA => B
X => Y, X => B, Z => ()
end) == ([StockAndFlow0([:A, :AA],[:B],[:A => :B, :AA => :B]), 
    StockAndFlow0([:X, :Z],[:Y, :B],[:X => :Y, :X => :B])])
```

---

# Stratification

## @stratify

Given three stockflows X, Z, Y and an expression block, create a new
stockflow representing the pullback of X -> Z and Y -> Z.

Equivalent to explicitly defining the StockFlow homomorphisms X -> Z and
Y -> Z and taking their pullback.

The particular morphisms can be inferred based on what each object in X and
Z maps to in Y.

For :stocks, :flows, :dynamic\_variables, :parameters and :sums in X and Y,
indicate what each element of X maps to on the left and Y on the right.  So,
if stocks x1, x2 in X and z in Z map to y in Y, you write 
:stocks
x1, x2 => y <= z

Prefix ~ to indicate a substring match and \_ to match everything else.

If there only exists one of a particular type of object which can be mapped
to, the maps don't need to be made explicit.  EG, if there only exists one
sum variable in Y, you don't need to have a :sums header.

If there exist multiple matches for an object (eg, you have A twice, or you
match A then have an \_), then only the first will be used.

Using stockflows with duplicate names could lead to unpredictable results
and is strongly recommended against.

```julia
@stratify WeightModel l_type ageWeightModel begin

    # Don't need stocks header, because there's only one
    :stocks
    _ => pop <= _

    :flows
    ~Death => f_death <= ~Death
    ~id => f_aging <= ~aging
    ~Becoming => f_fstOrder <= ~id
    _ => f_birth <= f_NB


    :dynamic_variables
    v_NewBorn => v_birth <= v_NB
    ~Death => v_death <= ~Death
    ~id  => v_aging <= v_agingCA, v_agingAS
    # Could match the right side with ~id
    v_BecomingOverWeight, v_BecomingObese => v_fstOrder <= v_idC, v_idA, v_idS

    :parameters
    μ => μ <= μ
    δw, δo => δ <= δC, δA, δS
    rw, ro => rFstOrder <= r
    rage => rage <= rageCA, rageAS

    # Similar to stocks, don't need a sums header.
    :sums
    N => N <= N

end 
```

---

Given n stockflows A\_1, ..., A\_{n-1}, Z and an expression block, create a 
new stockflow representing the pullback of A\_1 -> Z, ..., A\_{n-1} -> Z.

Use an ordered list of tuples to indicate the ith stockflow.  If there is a
single object in A\_i mapping to an object in Z, the tuple isn't necessary.
If 0 map to it, use an empty tuple (though, in that case, why does Z have
it at all?)

_ still acts as a default match, and ~ as a subtring match.  They each
act on their particular stockflow.

```julia
@n_stratify WeightModel ageWeightModel l_type begin

    # Once again, stocks header isn't necessary
    :stocks
    [_, _] => pop

    :flows
    [~Death, ~Death] => f_death
    [~id, ~aging] => f_aging 
    [~Becoming, ~id] => f_fstOrder
    [_, f_NB] => f_birth


    :dynamic_variables
    [v_NewBorn, v_NB] => v_birth
    [~Death, ~Death] => v_death
    [~id, (v_agingCA, v_agingAS)] => v_aging
    [(v_BecomingOverWeight, v_BecomingObese), (v_idC, v_idA, v_idS)] => v_fstOrder

    :parameters
    [μ, μ] => μ
    [(δw, δo), (δC, δA, δS)] => δ
    [(rw, ro), r] => rFstOrder
    [rage, (rageCA, rageAS)] => rage

    # Nor is sums header necessary
    :sums
    [N,N] => N
end
```

---

# Composition

## @compose

Given n stockflows A\_1, ..., A\_n and an expression block, return a new
stockflow such that specified shared stocks, sums and links between them
are treated as the same.

First line of the expression block must be the aliases used for the
corresponding stockflow in the block.  Every alias must be unique.

Use the same syntax for @foot to specify feet to be composed on.  Each line
takes the form X, Y ^ A => B, C => (), where the left side of ^ are
stockflows, the right are feet.

Cannot compose on an empty foot.  Cannot use a foot which has been used
on a previous line.  Changing the order of stock-sum links - eg, 
A => B, C => () and C => (), A => B - will be treated as different feet.
Each stockflow involved in a particular composition must have all stocks
and sums which are in the corresponding foot.

Composition with no stockflows given as argument returns an empty stockflow.

```julia
(@compose begin 
    ()
end) == (@compose begin end) == StockAndFlow()

sirv = @compose sir svi begin
  (sir, svi)
  sir, svi ^ S => N, I => N
end

XAY_model = @compose X SIS_A SIS_Y begin
    (X, A, Y)
    X, A, Y ^ X => N
    A, Y ^ () => NI
end

Diabetes_Model = @compose Model_Normoglycemic Model_Hyperglycemic Model_Norm_Hyper begin
    (Normo, Hyper, NH)

    Normo, NH ^ NormalWeight => N
    Normo, NH ^ OverWeight => N
    Normo, NH ^ Obese => N

    Hyper, NH ^ Prediabetic_U => N
    Hyper, NH ^ Prediabetic_D => N
end
```

---

# Rewrite

## @rewrite

Given a stockflow sf and an expression block, create a new stockflow with
edits made based on the block.  The expression block will be used to create
three stockflows L, I, R such that sf ⊇ L ⊇ I and R ⊇ I, then apply a
rewrite rule to replace L with R in sf.

Every object in sf must have a unique name.

Rewrite has 8 headers, :stocks, :flows, :sums, :dynamic\_variables, :redefs, :removes,
and :dyvar\_swaps

The first five are used to indicate if an instance of that type is being
added.  Use the same definition syntax as in @stock_and_flow.

:redefs is used to change the definition of a dynamic variable, flow or
sum.  Again, use the same syntax as in @stock_and_flow.

:removes indicates that an object should be deleted.  You just need the name of it in the original.

:dyvar\_swaps is used to swap all instances of an object inside dynamic variables
with another object.  Use the notation A => B to indicate every A in a
dynamic variable should now instead be B.

For this to work, there must exist homomorphisms I => L and I => R, and every name in the original stockflow must be unique.

```julia
aged_sir_rewritten = @rewrite aged_sir begin

    :redefs
    v_meanInfectiousContactsPerSv_cINC = cc_C * v_prevalencev_INC_post
    v_meanInfectiousContactsPerSv_cINA = cc_A * v_prevalencev_INA_post


    :parameters
    fcc
    fca
    fac
    faa

    :dynamic_variables

    v_CCContacts = fcc * v_prevalencev_INC
    v_CAContacts = fca * v_prevalencev_INA

    v_ACContacts = fac * v_prevalencev_INC
    v_AAContacts = faa * v_prevalencev_INA

    v_prevalencev_INC_post = v_CCContacts + v_CAContacts
    v_prevalencev_INA_post = v_ACContacts + v_AAContacts

end


Covid19_rewritten = @rewrite COVID19 begin
  :dyvar_swaps
  λ => v_NewIncidence₂
  rw_v => rw

  :removes
  λ
  rw_v
end

sirv_rewritten = @rewrite sirv begin
  :dyvar_swaps
  lambda => v_inf₂
  rdeath_svi => rdeath

  :removes
  lambda
  rdeath_svi
end
```

---

# Homomorphism

## @hom

Using two stockflows and a block with the same headers as @stratify, define an ACSetTransformation from the first stockflow to another.

Half of a stratify.  Rather than creating two homomorphisms and taking the pullback, creates one homomorphism.

```julia
@hom WeightModel l_type begin
  :stocks
  NormalWeight => pop
  OverWeight => pop
  Obese => pop

  :parameters
  μ => μ
  ~δ => δ
  rage => rage
  _ => rFstOrder

  :dynamic_variables
  v_NewBorn => v_birth
  ~Becoming => v_fstOrder
  ~Death => v_death
  _ => v_aging

  :flows
  f_NewBorn => f_birth
  ~Becoming => f_fstOrder
  ~Death => f_death
  ~id => f_aging
  _ => f_death

end
```
