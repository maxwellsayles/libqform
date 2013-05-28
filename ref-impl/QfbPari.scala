/**
 * This is NUCOMP from Pari.  The source is from basemath/Qfb.c
 * and has been ported to Scala.
 *
 * Specifically I was interested in prototyping different xgcd and pxgcd
 * functions and scala makes it easier to prototype than C.
 */

case class Qfb(a: BigInt, b: BigInt, c: BigInt)

def reduce(form: Qfb): Qfb = {
  form match {
    case Qfb(a, b, c) if (a == c && b < 0) => Qfb(a, -b, c)
    case Qfb(a, b, c) if (a > c) => reduce(Qfb(c, -b, a))
    case Qfb(a, b, c) if (b > a || b <= -a) =>
      val delta = b * b - 4 * a * c
      val a2 = 2 * a
      val bp = b % a2
      val bpp = if (bp > a) bp - a2 else if (b <= -a) bp + a2 else bp
      val cp = (bpp * bpp - delta) / (4 * a)
      reduce(Qfb(a, bpp, cp))
    case _ => form
  }
}

def xgcd(a: BigInt, b: BigInt): (BigInt, BigInt, BigInt) = {
  def f(u1: BigInt, u2: BigInt, u3: BigInt,
	v1: BigInt, v2: BigInt, v3: BigInt): (BigInt, BigInt, BigInt) = {
    if (v3 == 0) (u1, u2, u3)
    else {
      val q: BigInt = u3 / v3
      f(v1, v2, v3, u1 - q * v1, u2 - q * v2, u3 - q * v3)
    }
  }
  val (u1, u2, u3) = f(1, 0, a.abs, 0, 1, b.abs)
  val u1p = if (a < 0) -u1 else u1
  val u2p = if (b < 0) -u2 else u2
  assert(u1p * a + u2p * b == u3)
  assert(u3 > 0)
  (u1p, u2p, u3)
}

/*
def pxgcd(in_a: BigInt, in_b: BigInt, bound: BigInt):
	 (BigInt, BigInt, BigInt, BigInt, Int) = {
  assert(bound >= 0)
  var a: BigInt = in_a
  var b: BigInt = in_b
  var v: BigInt = 0
  var t: BigInt = 1
  var z: Int = 0
  while (b != 0 && b.abs > bound) {
    val q = a / b
    var r = a % b
    a = b
    b = r
    r = v - q * t
    v = t
    t = r
    z += 1
  }
  (a, b, v, t, z)
}
*/

def numBits(a: BigInt): Int = {
  assert(a != 0)
  a match {
    case _ if a < 0 => numBits(-a)
    case _ if a == 1 => 1
    case _ => 1 + numBits(a >> 1)
  }
}

/** Simplified left-to-right binary partial xgcd */
def pxgcd(a: BigInt, b: BigInt, bound: BigInt):
         (BigInt, BigInt, BigInt, BigInt, Int) = {
  assert(bound >= 0)
  assert(a.abs >= b.abs)
  def f(r1: BigInt, s1: Int,
	r0: BigInt, s0: Int,
	c1: BigInt, c0: BigInt,
	z: Int):
       (BigInt, BigInt, BigInt, BigInt, Int) = {
    if (r1 < 0) f(-r1, -s1, r0, s0, -c1, c0, z)
    else if (r0 < 0) f(r1, s1, -r0, -s0, c1, -c0, z)
    else if (r1 < r0) f(r0, s0, r1, s1, c0, c1, z + 1)
    else if (r0 == 0 || r0 <= bound) {
      (s1 * r1, s0 * r0, s1 * c1, s0 * c0, z)
    } else {
      val k = numBits(r1) - numBits(r0)
      val rp = r1 - (r0 << k)
      val cp = c1 - (c0 << k)
      f(rp, s1, r0, s0, cp, c0, z)
    }
  }
  f(a, 1, b, 1, 0, 1, 0)
}

def nucomp(form1: Qfb, form2: Qfb): Qfb = {
  var Qfb(a1, b1, c1) = form1
  var Qfb(a2, b2, c2) = form2
  assert(a1 > 0)
  assert(a2 > 0)
  var delta: BigInt = b1 * b1 - 4 * a1 * c1
  var delta4: BigInt = Math.pow(delta.doubleValue.abs, 1.0/4.0).toInt
  if (a1.abs < a2.abs) {
    var t = a1
    a1 = a2
    a2 = t
    t = b1
    b1 = b2
    b2 = t
    t = c1
    c1 = c2
    c2 = t
  }

  var s = (b1 + b2) / 2
  val n = b2 - s
  var (u, v, d) = xgcd(a2, a1)
  var a: BigInt = 0
  var b: BigInt = 0
  var d1: BigInt = 0
  var p1: BigInt = 0
  var p2: BigInt = 0
  var a3: BigInt = 0
  var b3: BigInt = 0
  var c3: BigInt = 0
  var g: BigInt = 0
  if (d == 1) {
    a = -u * n
    d1 = d
  } else {
    val (u1, v1, d1p) = xgcd(s, d)
    d1 = d1p
    if (d1 != 1) {
      a1 /= d1
      a2 /= d1
      s /= d1
      d /= d1
    }
    p1 = c1 % d
    p2 = c2 % d
    val l = (-u1 * (u * p1 + v * p2)) % d
    a = (l * a1 / d) - (u * n / d)
  }

  a %= a1
  d = a1
  var (dp, v3, vp, v2, z) = pxgcd(d, a, delta4)
  d = dp
  v = vp
  if (z == 0) {
    // Normal composition.
    g = (v3 * s + c2) / d
    b = a2
    v2 = d1
    a3 = d * b
  } else {
    // NUCOMP.
    if ((z & 1) == 1) {
      v3 = -v3
      v2 = -v2
    }
    b = (a2 * d + n * v) / a1
    val e = (s * d + c2 * v) / a1
    val q3 = e * v2
    val q4 = q3 - s
    b2 = q3 + q4
    g = q4 / v
    assert(d1 > 0)
    if (d1 != 1) {
      v2 *= d1
      v  *= d1
      b2 *= d1
    }
    a3 = d * b + e * v
  }
  val q1 = b * v3
  val q2 = q1 + n;
  b3 = b2 + (if (z > 0) q1 + q2 else 2 * q1)
  c3 = (v3 * q2 / d) + (g * v2)
  reduce(Qfb(a3, b3, c3))
}

val q = Qfb(11, 4, 5772715433271L)
//val q = Qfb(11, 9, BigInt("48347978850011094692042314858"))
println("1 -> " + q)
var qp = q
var i = 2
while (qp.a != 1) {
  qp = nucomp(qp, q)
  println("%d -> %s".format(i, qp))
  i += 1
}

