NAME                ILLCOND1
ROWS
 N  COST
 G  C1
 G  C2
 G  C3
 G  C4
COLUMNS
    X1              COST              1       C1                105
    X1              C2              100       C3                3
    X1              C4                1
    X2              COST              1       C1                100
    X2              C2              105       C3                1
    X2              C4                3
RHS
    rhs             C1              150       C2                150
    rhs             C3                3       C4                3
BOUNDS
 UI BOUND           X1                5
 UI BOUND           X2                5
ENDATA