
#= ---------------------------------------------------------------------
Metabolites
ATP #1, l-citrulline #2, l-asparate #3, AMP #4,
diphosphate #5, (Nomega-L-arginino)succinate #6, fumarate #7,
l-arginine #8, H2O #9, ornithine #10, urea #11,
carbamoyl-phosphate #12, phosphate #13, NADPH #14, H #15, O2 #16,
Nitric oxide #17, NADP #18
Reactions
v1, v2, v3, v4, v5fwd, v5frev
b1- carbomyl-phosphate in, b2- aspartate in, b3- fumarate out,
b4- urea out, b5- ATP in, b6- AMP out, b7- diphosphate out,
b8- water in, b9- phosphate out, b10- NADPH in, b11- H+ in,
b12- O2 in, b13- NO out, b14- NADP+ out, b15- water out
Atom array
C, H, N, O, P, S
----------------------------------------------------------------------=#

using CSV
using LinearAlgebra
using DataFrames

stoichiometric_matrix= DataFrame(CSV.File("Network.csv", header=false))
stoichiometric_matrix= convert(Matrix, stoichiometric_matrix)
elemental_matrix= DataFrame(CSV.File("Atom.csv", header=false))
elemental_matrix= convert(Matrix, elemental_matrix)

check_atom_balance = transpose(elemental_matrix)*stoichiometric_matrix

nrows, ncols = size(check_atom_balance)
atom_balance_array = zeros(nrows)

for row in collect(1:nrows)
    atom_balance_array[row] = sum((check_atom_balance[row,:]))
end

return atom_balance_array
#Atoms are balanced if all elements are equal to zero.
#Nitrogen balance is negative
