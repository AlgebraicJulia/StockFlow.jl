module Rewrite


export sfrewrite, @rewrite

struct SFType
  type::Symbol
  name::String
  index::Int
end



using ...StockFlow
using ..StockFlow.Syntax
import ..StockFlow.Syntax: parse_dyvar, parse_flow, parse_sum, parse_stock, 
parse_param, sf_to_block, StockAndFlowBlock, stock_and_flow_syntax_to_arguments
using MLStyle


function reset_positions(sfold, sfnew)
  # 1: select all links to dv which have a position:
  oldlv = Dict((i, j) => k for (i, j, k) ∈ zip(get_lvs(sfold), get_lvv(sfold), get_lvsposition(sfold)))
  oldlsv = Dict((i, j) => k for (i, j, k) ∈ zip(get_lsvsv(sfold), get_lsvv(sfold), get_lsvsvposition(sfold)))
  oldlvv = Dict((i, j) => k for (i, j, k) ∈ zip(get_lvsrc(sfold), get_lvtgt(sfold), get_lvsrcposition(sfold)))
  oldlpv = Dict((i, j) => k for (i, j, k) ∈ zip(get_lpvp(sfold), get_lpvv(sfold), get_lpvpposition(sfold)))

  oldstocks = snames(sfold)
  olddyvars = vnames(sfold)
  oldparams = pnames(sfold)
  oldsums = svnames(sfold)
  
  newlv = zip(get_lvs(sfnew), get_lvv(sfnew), get_lvsposition(sfnew))
  newlsv = zip(get_lsvsv(sfnew), get_lsvv(sfnew), get_lsvsvposition(sfnew))
  newlvv = zip(get_lvsrc(sfnew), get_lvtgt(sfnew), get_lvsrcposition(sfnew))
  newlpv = zip(get_lpvp(sfnew), get_lpvv(sfnew), get_lpvpposition(sfnew))

  newstocks = snames(sfnew)
  newdyvars = vnames(sfnew)
  newparams = pnames(sfnew)
  newsums = svnames(sfnew)
  
  #2: iterate over all new values, find what they link to in new, use that to find the corresponding symbol in old

  restored_lvsposition = zeros(Int, nlv(sfnew))
  restored_lsvsvposition = zeros(Int, nlsv(sfnew))
  restored_lvsrcposition = zeros(Int, nlvv(sfnew))
  restored_lpvpposition = zeros(Int, nlpv(sfnew))
  
  for (lv, (lvs, lvv, lvspos)) ∈ enumerate(newlv)
      old_pos = oldlv[(lvs, lvv)]
      restored_lvsposition[lv] = old_pos
  end

      
  for (lsv, (lsvsv, lsvv, lsvsvpos)) ∈ enumerate(newlsv)
      old_pos = oldlsv[(lsvsv, lsvv)]
      restored_lsvsvposition[lsv] = old_pos
  end

  
  for (lvv, (lvsrc, lvtgt, lvsrcpos)) ∈ enumerate(newlvv)
      old_pos = oldlvv[(lvsrc, lvtgt)]
      restored_lvsrcposition[lvv] = old_pos
  end

  for (lpv, (lpvp, lpvv, lpvppos)) ∈ enumerate(newlpv)
      old_pos = oldlpv[(lpvp, lpvv)]
      restored_lpvpposition[lpv] = old_pos
  end

  set_subpart!(sfnew, :lvsposition, restored_lvsposition)
  set_subpart!(sfnew, :lsvsvposition, restored_lsvsvposition)
  set_subpart!(sfnew, :lvsrcposition, restored_lvsrcposition)
  set_subpart!(sfnew, :lpvpposition, restored_lpvpposition)

  return sfnew
end


function sfrewrite(sf::K, block::Expr) where {K < AbstractStockAndFlowF}
  Base.remove_linenums!(block)


  

  sf_snames = invert_vector(snames(sf))
  sf_svnames = invert_vector(svnames(sf))
  sf_vnames = invert_vector(vnames(sf))
  sf_fnames = invert_vector(fnames(sf))
  sf_pnames = invert_vector(pnames(sf))



  # # May be unnecessary.
  name_vector = [snames(sf) ; svnames(sf) ; vnames(sf) ; fnames(sf) ; pnames(sf)]
  @assert allunique(name_vector) "Not all names are unique!  $(filter(x -> count(y -> y == x, name_vector) >= 2, name_vector))"



  sf_sdict::Dict{Int::SFType} = Dict{Int::SFType}(i => SFType(:S, sf_snames[i], i) for i ∈ eachindex(sf_snames))
  sf_svdict::Dict{Int::SFType} = Dict{Int::SFType}(i => SFType(:SV, sf_svnames[i], i) for i ∈ eachindex(sf_svnames))
  sf_vdict::Dict{Int::SFType} = Dict{Int::SFType}(i => SFType(:V, sf_vnames[i], i) for i ∈ eachindex(sf_vnames))
  sf_fdict::Dict{Int::SFType} = Dict{Int::SFType}(i => SFType(:F, sf_fnames[i], i) for i ∈ eachindex(sf_fnames))
  sf_pdict::Dict{Int::SFType} = Dict{Int::SFType}(i => SFType(:P, sf_pnames[i], i) for i ∈ eachindex(sf_pnames))





  sf_allnames::Dict{Int::SFType} = merge(sf_sdict, sf_svdict, sf_vdict, sf_fdict, sf_pdict)


  sf_block::StockAndFlowBlock = sf_to_block(sf)

  


  # So we have all the names and the indices now.
  # I: Intermediate.  Common parts.
  # L: Old.  Parts that including parts which are changing.
  # R: New.

  # +: Add to R.
  # -: Add to L.
  # :swap x => y: Add x to L.  Add everything which contains x to L.
  #   Add everything which contains x, with x replaced with y, to R.
  #   Add everything which contains x, minus the x itself, to I.

  # Ignore the following: 
  # :swap x := y: Add x to L.  Add everything which contains x to L.
  #   Add the new definition of x to R.
  #   Add difference of R - L to I.

  # The former allows to implictly update all instances

  # Latter swap allows for redefinition, former just replaces.

  # We just need the former, though.




end





macro rewrite(sf, block)
  escaped_block = Expr(:quote, block)
  sf = esc(sf)
  quote
    sfrewrite($sf, $escaped_block)
  end
end





end