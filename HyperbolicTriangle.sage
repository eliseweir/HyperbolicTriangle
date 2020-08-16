"""
Creates nxn matrices a, b, c satisfying a^p = b^q = c^r = abc = id.

Produces a representation of a T(p,q,r) hyperbolic triangle group in PSL(n,R) by composing the holonomy 
representation T(p,q,r) -> PSL(2,R) with the unique irreducible representation PSL(2,R) -> PSL(n,R).
"""

import sys

class HyperTriAngle:
    """
    Represents an angle in a hyperbolic triangle and the side opposite that angle.
    """
    def __init__(self,order):
        self.order = order
        self.angle = pi/order

def homogeneous(n,x,y):
    """
    Create a homogeneous degree n - 1 polynomial, e.g., f0 x^2 + f1 x y + f2 y^2 for n = 3
    """
    poly = 0
    for k in range(0,n):
        fsubk = var('f'+str(k))
        poly += fsubk * x**(n - 1 - k) * y**k
    return poly

def irreducible(n,rep):
    """
    Find a representation for rep in SL(n,R) under the unique irreducible representation SL(2,R) -> SL(n,R)
    """
    # use the 2x2 representation rep to create a variable substitution
    x0 = vector([x, y])* (vector([1, 0]) * rep)
    y0 = vector([x, y])* (vector([0, 1]) * rep)
    poly = homogeneous(n,x0,y0).expand()

    # extract coefficients from the polynomial using derivatives
    repn = Matrix(SR,n,n)
    for i in range(0,n):
        print('.',end="")
        for j in range(0,n):
            repn[i,j] = derivative(derivative(derivative(poly,var('f'+str(j))),x,n - i - 1),y,i)/(factorial(n-i-1)*factorial(i))

    return repn

def prettyprint(mat,file):
    """
    Print matrices to file in a format that can be copied directly into Sagemath (or similar)
    """
    file.write('[')
    # print each row
    for i in range(0,n-1):
        file.write('[')
        # print each column, separated by commas
        for j in range (0,n-1):
            file.write(str(mat[i][j]) + ', ')
        # the last entry in each row shouldn't be followed by a comma
        file.write(str(mat[i][n-1]))
        file.write('], ')
    # the last row shouldn't end with a comma
    file.write('[')
    for k in range (0,n-1):
        file.write(str(mat[n-1][k]) + ', ')
    file.write(str(mat[n-1][n-1]))
    file.write(']')
    file.write(']\n')

# accept angles pi/p, pi/q, pi/r for a hyperbolic triangle
angleflag = 1
while angleflag:
    pqr_raw = input("Input your T(p,q,r) triangle as p,q,r: ")
    pqr = pqr_raw.split(',')
    if len(pqr) != 3:
        print("You must enter 3 angles.\n")
    else:
        try:
            p = int(pqr[0])
            q = int(pqr[1])
            r = int(pqr[2])
        except:
            print("\nError: p,q,r must be integers.")
            sys.exit()
        # check to make sure the triangle is actually hyperbolic
        if p < 2 or q < 2 or r < 2:
            print("p,q,r must be at least 2 so that angles are at most pi/2.\n")
        elif 1/p + 1/q + 1/r >= 1:
            print("Your T(p,q,r) triangle must be hyperbolic. You entered a", end = " ")
            if 1/p + 1/q + 1/r == 1:
                print("Euclidean triangle.\n")
            elif 1/p + 1/q + 1/r >= 1:
                print("spherical triangle.\n")
        else:
            angleflag = 0

# accept the size of matrices to output
dimensionflag = 1
while dimensionflag:
    try:
        n = int(input("Enter the degree n of the nxn matrices you want as outputs: "))
    except:
        print("\nError: For nxn matrices, n must be an integer.")
        sys.exit()
    # check to make sure the size of matrices requested is at least 2
    if n < 2:
        print("nxn matrices must have size at least n = 2.\n")
    # warn about computing time for large matrices
    elif n > 7:
        yesno = input("Matrices larger than 7x7 may take a long time. Are you sure? Y/N: ")
        if yesno.upper() == 'Y' or yesno.upper() == 'YES':
            dimensionflag = 0
        else:
            print("Exiting...")
            sys.exit()
    else:
        dimensionflag = 0

a = HyperTriAngle(p)
b = HyperTriAngle(q)
c = HyperTriAngle(r)

# matrices for reflection over sides opposite b and c
b.reflection = (Matrix([[1, cos(a.angle)], [0, 1]]) * Matrix([[0, 1], [1, 0]])) * Matrix([[1, -cos(a.angle)], [0, 1]])
c.reflection = Matrix([[-1, 0], [0, 1]])

# use hyperbolic geometry to find the length of the side opposite c and use it to find the coordinate
# along the imaginary axis
c.length = (cos(a.angle)*cos(b.angle) + cos(c.angle))/(sin(a.angle) * sin(b.angle))


# move along the imaginary axis from i by the length of the side opposite c
var('y')
sols = solve((c.length == 1 + ((0 - 0)**2 + (y - sin(a.angle))^2)/(2 * sin(a.angle) * y)),y)
c.axis = sols[1].rhs()

# find the radius and horizontal displacement from the origin of a circle making angle b with imaginary axis
b.radius = c.axis/sin(b.angle)
b.disp = b.radius * cos(b.angle)

# use the radius and displacement to translate and scale a reflection over the unit circle
translA = Matrix([[1, b.disp], [0, 1]])
scaleA = Matrix([[1/sqrt(b.radius),0],[0,sqrt(b.radius)]])
a.reflection = (~translA * ~scaleA * Matrix([[0, 1], [1, 0]]) * scaleA * translA).simplify_trig()

# compose reflections over sides opposite a, b, and c to get 2x2 representations
a.rep2x2 = b.reflection * c.reflection
b.rep2x2 = c.reflection * a.reflection
c.rep2x2 = a.reflection * b.reflection

# use a variable substitution in homogeneous polynomials to get nxn representations from 2x2 matrices
print("Computing a.",end='')
a.repnxn = (irreducible(n,a.rep2x2)).simplify_full()
print(" done.\nComputing b.",end='')
b.repnxn = (irreducible(n,b.rep2x2)).simplify_full()
print(" done.\nComputing c.",end='')
c.repnxn = (irreducible(n,c.rep2x2)).simplify_full()
print(" done.")

# save the results, both exact and approximate, to hyperbolictriangle_output.txt
with open('hyperbolictriangle_output.txt','w') as output_file:
    output_file.write('Exact representations of T(p,q,r) = <a,b,c| a^p = b^q = c^r = abc = 1>:\n')
    output_file.write('\na = ')
    prettyprint(a.repnxn,output_file)
    output_file.write('\nb = ')
    prettyprint(b.repnxn,output_file)
    output_file.write('\nc = ')
    prettyprint(c.repnxn,output_file)
    output_file.write('\nNumerical approximations:\n')
    output_file.write('\na = ')
    prettyprint((a.repnxn.n()),output_file)
    output_file.write('\nb = ')
    prettyprint((b.repnxn.n()),output_file)
    output_file.write('\nc = ')
    prettyprint((c.repnxn.n()),output_file)

print("Output written to hyperbolictriangle_output.txt in current directory.")
