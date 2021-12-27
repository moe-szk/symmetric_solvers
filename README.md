# symmetric_solvers

Liner equation solvers for sparse symmetric coefficient matrix
 
### Solvers without preconditioning 
- CG
- MINRES
- CR
- MR
- GD
- SOR


### Solvers with preconditioning (incomplete Cholesky decomposition)
- Preconditioned CG (=ICCG)
- Preconditioned MINRES
- Preconditioned CR
- Preconditioned MR
- Preconditioned GD

### CG with variable preconditioning
- ICCG
- MINRES
- CR

### Tools 
- m_icdc1d

incomplete Cholesky decomposition

- m_conv_2_iccg

compression

- m_prod1d

$$ A {\bf x} $$

- m_icsl1d

$$ A^{-1} {\bf x} $$
