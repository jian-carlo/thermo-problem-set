import sympy as sp
from typing import Union

def newtons_method_debug(function:sp.Expr, variable:sp.Symbol, guess: Union[float, int], tolerance=1e-7, max_iterations=1000):
    """
    Newton's Method to numerically solve for the root of an algebraic function.
    `newtons_method` is a version of this function without print statements.

    Parameters:
    function (Expr): The SymPy expression representing the function.
    variable (Symbol): The independent variable of the function.
    guess (float or int): Initial guess for the root.
    tolerance (float, optional): The accuracy of the solution (default is 1e-7).
    max_iterations (int, optional): Maximum number of iterations (default is 1000).

    Returns:
    float: The approximated root of the function.

    Raises:
    ValueError: If the derivative is near zero or if the method fails to converge.
    """
    x = variable

    f = function
    df = sp.diff(function)

    f_ = sp.lambdify(x, f)
    df_ = sp.lambdify(x, df)

    xn = guess
    for n in range(0, max_iterations):

        fxn = f_(xn)

        if abs(fxn) < tolerance:
            print(f"Found solution after {n} iterations.")
            print(f"root = {xn}")
            return xn

        dfxn = df_(xn)

        if dfxn == 0:
            print("Zerp derivative. No solution found.")
            return None

        xn = xn - fxn/dfxn
    print("Exceeded maximum iterations. No solution found.")
    return None 

def newtons_method(function:sp.Expr, variable:sp.Symbol, guess: Union[float, int], tolerance=1e-7, max_iterations=1000):
    """
    Newton's Method to numerically solve for the root of an algebraic function.

    Parameters:
    function (Expr): The SymPy expression representing the function.
    variable (Symbol): The independent variable of the function.
    guess (float or int): Initial guess for the root.
    tolerance (float, optional): The accuracy of the solution (default is 1e-7).
    max_iterations (int, optional): Maximum number of iterations (default is 1000).

    Returns:
    float: The approximated root of the function.

    Raises:
    ValueError: If the derivative is near zero or if the method fails to converge.
    """
    x = variable

    f = function
    df = sp.diff(function)

    f_ = sp.lambdify(x, f)
    df_ = sp.lambdify(x, df)

    xn = guess
    for n in range(0, max_iterations):

        fxn = f_(xn)

        if abs(fxn) < tolerance:
            return xn

        dfxn = df_(xn)

        if dfxn == 0:
            return None

        xn = xn - fxn/dfxn
    return None 
