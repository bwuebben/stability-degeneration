"""Independent r02 certificate using integer intervals and fixed-frame derivatives.

No numerical library, eigensolver, moving singular frame, or finite difference
is used.  Endpoints are integer multiples of 10**-60.  Every arithmetic
operation rounds outward using integer division; square roots use math.isqrt.
Exponential Frechet derivatives are summed in their original fixed frames,
with norm bounds for all omitted terms.  Run with the Python standard library.
"""
from fractions import Fraction
from math import factorial, isqrt
from pathlib import Path
import hashlib
import json

SCALE = 10**60


class I:
    def __init__(self, value=0, *, endpoints=None):
        if endpoints is not None:
            self.lo, self.hi = endpoints
        elif isinstance(value, I):
            self.lo, self.hi = value.lo, value.hi
        else:
            f = Fraction(value) * SCALE
            self.lo = f.numerator // f.denominator
            self.hi = -((-f.numerator) // f.denominator)
        assert self.lo <= self.hi

    def __add__(self, other):
        other = I(other)
        return I(endpoints=(self.lo + other.lo, self.hi + other.hi))

    __radd__ = __add__

    def __neg__(self):
        return I(endpoints=(-self.hi, -self.lo))

    def __sub__(self, other):
        return self + -I(other)

    def __rsub__(self, other):
        return I(other) + -self

    def __mul__(self, other):
        other = I(other)
        values = [x*y for x in (self.lo, self.hi)
                  for y in (other.lo, other.hi)]
        return I(endpoints=(min(values)//SCALE, -((-max(values))//SCALE)))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = I(other)
        assert other.lo * other.hi > 0
        ratios = [Fraction(x*SCALE, y) for x in (self.lo, self.hi)
                  for y in (other.lo, other.hi)]
        low, high = min(ratios), max(ratios)
        return I(endpoints=(low.numerator//low.denominator,
                            -((-high.numerator)//high.denominator)))

    def sqrt(self):
        assert self.lo >= 0
        lower = isqrt(self.lo*SCALE)
        upper = isqrt(self.hi*SCALE)
        upper += upper*upper < self.hi*SCALE
        return I(endpoints=(lower, upper))

    def inflate(self, radius):
        radius = I(radius)
        assert radius.lo >= 0
        return I(endpoints=(self.lo-radius.hi, self.hi+radius.hi))

    def absolute_upper(self):
        return Fraction(max(abs(self.lo), abs(self.hi)), SCALE)

    def contains_zero(self):
        return self.lo <= 0 <= self.hi

    def __str__(self):
        def decimal(k):
            sign = '-' if k < 0 else ''
            integer, fraction = divmod(abs(k), SCALE)
            return f'{sign}{integer}.{fraction:060d}'
        return f'[{decimal(self.lo)}, {decimal(self.hi)}]'


def zeros(rows, cols=None):
    return [[I(0) for _ in range(cols or rows)] for _ in range(rows)]


def matrix(data):
    return [[I(x) for x in row] for row in data]


def eye(size):
    return matrix([[int(i == j) for j in range(size)] for i in range(size)])


def add(a, b):
    return [[x+y for x,y in zip(ar, br)] for ar,br in zip(a,b)]


def scale(a, scalar):
    return [[x*scalar for x in row] for row in a]


def mul(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a):
    return [list(row) for row in zip(*a)]


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def determinant(a):
    return a[0][0]*a[1][1]-a[0][1]*a[1][0]


def negative_projector(a):
    t = trace(a)
    gap = (t*t-4*determinant(a)).sqrt()
    return [[(I(int(i == j))-(2*a[i][j]-int(i == j)*t)/gap)/2
             for j in range(2)] for i in range(2)]


def exp_negative_rational(x, degree):
    assert x >= 0 and degree % 2 == 0
    even = sum(((-x)**k/Fraction(factorial(k)) for k in range(degree+1)), Fraction(0))
    odd = even + (-x)**(degree+1)/Fraction(factorial(degree+1))
    assert odd < even
    # The integral Taylor remainder has the sign of the next derivative.
    return I(endpoints=(I(odd).lo, I(even).hi))


def infinity_norm_bound(a):
    return max(sum(x.absolute_upper() for x in row) for row in a)


def exponential_jet(n, dn, degree=180):
    """Enclose exp(n), Dexp_n(dn), with explicit positive norm tails."""
    # Partial terms P_k=N^k/k!, Q_k=D(N^k/k!)[dN].
    p, q = eye(len(n)), zeros(len(n))
    ep, eq = p, q
    for k in range(1, degree+1):
        q = scale(add(mul(q,n), mul(p,dn)), Fraction(1,k))
        p = scale(mul(p,n), Fraction(1,k))
        ep, eq = add(ep,p), add(eq,q)
    m, ell = infinity_norm_bound(n), infinity_norm_bound(dn)
    assert m < degree+1
    # Geometric domination of successive terms, using submultiplicativity.
    tail = m**(degree+1)/factorial(degree+1)/(1-m/Fraction(degree+2))
    dtail = ell*m**degree/factorial(degree)/(1-m/Fraction(degree+1))
    ep = [[x.inflate(tail) for x in row] for row in ep]
    eq = [[x.inflate(dtail) for x in row] for row in eq]
    return ep,eq,dict(norm_bound=str(m),derivative_norm_bound=str(ell),
                      exponential_tail=str(I(tail)),derivative_tail=str(I(dtail)))


def put(target, row, col, block):
    for i,line in enumerate(block):
        for j,x in enumerate(line):
            target[row+i][col+j] = x


def block(a,row,col,rows=2,cols=2):
    return [line[col:col+cols] for line in a[row:row+rows]]


def show(a):
    return [[str(x) for x in row] for row in a] if isinstance(a,list) else str(a)


def main():
    b1,b2 = I(Fraction(1,25)),I(4)
    p = exp_negative_rational(Fraction(1,25),60)
    q = exp_negative_rational(Fraction(4),200)
    r1,r2 = (1-p)/b1,(1-q)/b2
    f1,f2 = r1-p,r2-q
    ell = (p-q)/(b2-b1)
    sigma = I(Fraction(2,5))*(r1-r2)/(b2-b1)
    h = matrix([[2*p,2*ell],[2*ell,5*q]])
    g = matrix([[f1-p,2*sigma],[2*sigma,4*f2-q]])
    w = 4-r1-4*r2
    assert determinant(h).hi < 0 and determinant(g).hi < 0 and w.lo > 0
    pv,pu = negative_projector(h),negative_projector(g)
    d,a = matrix([[Fraction(1,5),0],[0,2]]),matrix([[1,2]])
    dd = scale(add(mul(d,pv),mul(pu,d)),Fraction(1,2))
    da = scale(mul(a,pv),Fraction(1,2))
    dt,at = transpose(d),transpose(a)
    b,c = mul(dt,d),mul(at,a)
    db = add(mul(transpose(dd),d),mul(dt,dd))
    dl = add(mul(dd,dt),mul(d,transpose(dd)))
    dc = add(mul(transpose(da),a),mul(at,da))
    e = mul(mul(d,c),dt)
    de = add(add(mul(mul(dd,c),dt),mul(mul(d,dc),dt)),mul(mul(d,c),transpose(dd)))

    nv,dnv = zeros(4),zeros(4)
    put(nv,0,0,scale(b,-1)); put(nv,2,2,scale(b,-1)); put(nv,0,2,c)
    put(dnv,0,0,scale(db,-1)); put(dnv,2,2,scale(db,-1)); put(dnv,0,2,dc)
    ev,dev,bv = exponential_jet(nv,dnv)
    dh = add(block(dev,0,0),block(dev,0,2))
    hexp = add(block(ev,0,0),block(ev,0,2))

    nu,dnu = zeros(6),zeros(6)
    put(nu,0,0,scale(b,-1)); put(nu,4,4,scale(b,-1))
    put(nu,0,2,e); put(nu,2,4,eye(2))
    put(dnu,0,0,scale(dl,-1)); put(dnu,4,4,scale(dl,-1)); put(dnu,0,2,de)
    eu,deu,bu = exponential_jet(nu,dnu)
    dg = add(block(deu,0,4),scale(block(deu,0,0),-1))
    gexp = add(block(eu,0,4),scale(block(eu,0,0),-1))

    nr,dnr = zeros(4),zeros(4)
    put(nr,0,0,scale(b,-1)); put(nr,0,2,eye(2)); put(dnr,0,0,scale(db,-1))
    er,der,br = exponential_jet(nr,dnr)
    rb,drb = block(er,0,2),block(der,0,2)
    dw = -(mul(mul(da,rb),at)[0][0]+mul(mul(a,drb),at)[0][0]
           +mul(mul(a,rb),transpose(da))[0][0])
    trace_derivative = trace(dh)+trace(dg)+dw
    assert trace_derivative.contains_zero()
    for closed,series in [(h,hexp),(g,gexp)]:
        for i in range(2):
            for j in range(2):
                assert (closed[i][j]-series[i][j]).contains_zero()
    reduced = -2*(trace(mul(pv,dh))+trace(mul(pu,dg)))
    full = trace(mul(add(eye(2),scale(pv,-2)),dh))
    full += trace(mul(add(eye(2),scale(pu,-2)),dg))+dw
    for result in (reduced,full):
        assert result.lo > I('0.23957510898943').hi
        assert result.hi < I('0.23957510898945').lo
    result = {key:show(value) for key,value in dict(p=p,q=q,detH=determinant(h),
        detG=determinant(g),w=w,Pv=pv,Pu=pu,dD=dd,dA=da,dH_fixed=dh,dG_fixed=dg,
        dw=dw,trace_derivative=trace_derivative,derivative_over_pi=reduced,
        full_derivative_over_pi=full).items()}
    result['tails'] = dict(V=bv,U=bu,R=br)
    result['decimal_endpoint_scale'] = 60
    result['taylor_degree'] = 180
    result['certificate'] = '0.23957510898943 < limiting derivative / pi < 0.23957510898945'
    result['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    output = Path(__file__).parent/'results'/'mass-increase-review.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ('derivative_over_pi',
                    'full_derivative_over_pi','trace_derivative','certificate')},indent=2))


if __name__ == '__main__':
    main()
