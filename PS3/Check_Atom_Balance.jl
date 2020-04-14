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
