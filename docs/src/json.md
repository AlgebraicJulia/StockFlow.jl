# # Serializing Stock and Flow Models.


```@example JSON
using StockFlow
using Catlab
using Catlab.CategoricalAlgebra
using JSON3
import StockFlow: StockAndFlowpUntyped, vectorify, state_dict, add_stocks!, add_flow!, add_links!
```

In order to serialize our models to JSON, we need to encode our functions as strings instead of julia functions. You can always recover the original model with `Base.Meta.parse` and `eval`.

```@example JSON
StockAndFlowSymbolic = StockAndFlowpUntyped{Symbol,String} 

# This is our constructor specialized to String formulas.
NewStockAndFlowSymbolic(s,f) = begin
	d = StockAndFlowSymbolic()

    s = vectorify(s)
    add_stocks!(d,length(s),sname=s)

    s_idx = state_dict(s)
    f = vectorify(f)
    for (i, ((fattr,uds),ls)) in enumerate(f)
      fn = first(fattr)
      ff = last(fattr)
    	sui = s_idx[first(uds)]
    	sdi = s_idx[last(uds)]
    	ls = vectorify(ls)
    	add_flow!(d,sui,sdi,fname=fn,ϕf=ff)
    	add_links!(d,map(x->s_idx[x],ls),repeat([i], length(ls)), length(ls))
    end
    d
end
```

## Examples

Our first example model is SIR.

```@example JSON
# define functions ϕ of flows in the SIR model
fNewIncidence = "p.cβ*u.S*u.I/p.N"
fNewRecovery = "u.I/p.tr"

# StockAndFlowp(stocks,
#               (flow=>function, upstream=>downstream) => stocks linked)
sir = NewStockAndFlowSymbolic((:S, :I, :R), 
   ((:NewIncidence=>fNewIncidence, :S=>:I)=>(:S,:I),
    (:NewRecovery=>fNewRecovery, :I=>:R)=>:I))
```

```@example JSON
JSON3.pretty(generate_json_acset(sir))
```

We can include some of the premade models, which are found in the `PremadeModels` module. They are all full-bore Stock and Flows.

```@example JSON
m = StockFlow.PremadeModels.seir_model
```
```@example JSON
JSON3.pretty(generate_json_acset(m))
```

## Schemas

The restricted class of primitive models.

```@example JSON
JSON3.pretty(generate_json_acset_schema(TheoryStockAndFlowp))
```

And for the generic schema

```@example JSON
JSON3.pretty(generate_json_acset_schema(TheoryStockAndFlow))
```
