Ver 6.1:

made several changes to PTX_Tx:
1) for a given symbol, we walk from optimal b at Hamming range r from its predecessor to another  optimal b ar Hamming range r
2) corrected to possible PTS vector to such that really contains all the possible combinations, including identical phases, such as [1,1] or [j, j]
3) improved efficienc of code to skip the iterations where the assignment of the r terms to b does not change it

4) improved efficiency of PAPR_calc

5) changed the chunk of the signal on which the CCDF curve analyzes to not contain only data symbols




