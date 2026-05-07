/**
 * Exact Fraction arithmetic for the browser.
 * Mirrors Python's fractions.Fraction.
 */
class Frac {
  constructor(num = 0, den = 1) {
    if (den === 0) throw new Error("Zero denominator");
    if (typeof num === "string") {
      if (num.includes("/")) {
        const [a, b] = num.split("/").map(Number);
        num = a; den = b;
      } else {
        num = Number(num); den = 1;
      }
    }
    if (den < 0) { num = -num; den = -den; }
    const g = Frac.gcd(Math.abs(num), Math.abs(den));
    this.n = num / g;
    this.d = den / g;
  }

  static gcd(a, b) { a = Math.abs(a); b = Math.abs(b); while (b) { [a, b] = [b, a % b]; } return a || 1; }
  static from(v) { return v instanceof Frac ? v : new Frac(typeof v === "string" ? v : Number(v)); }
  static ZERO = new Frac(0);
  static ONE = new Frac(1);

  add(o) { o = Frac.from(o); return new Frac(this.n * o.d + o.n * this.d, this.d * o.d); }
  sub(o) { o = Frac.from(o); return new Frac(this.n * o.d - o.n * this.d, this.d * o.d); }
  mul(o) { o = Frac.from(o); return new Frac(this.n * o.n, this.d * o.d); }
  div(o) { o = Frac.from(o); if (o.n === 0) throw new Error("Division by zero"); return new Frac(this.n * o.d, this.d * o.n); }
  neg() { return new Frac(-this.n, this.d); }

  eq(o) { o = Frac.from(o); return this.n === o.n && this.d === o.d; }
  lt(o) { o = Frac.from(o); return this.n * o.d < o.n * this.d; }
  gt(o) { o = Frac.from(o); return this.n * o.d > o.n * this.d; }
  lte(o) { return this.lt(o) || this.eq(o); }
  gte(o) { return this.gt(o) || this.eq(o); }
  isZero() { return this.n === 0; }
  isNeg() { return this.n < 0; }
  isInt() { return this.d === 1; }

  floor() { return Math.floor(this.n / this.d); }
  fracPart() { return this.sub(new Frac(this.floor())); }

  toString() { return this.d === 1 ? `${this.n}` : `${this.n}/${this.d}`; }
  toHTML() {
    if (this.d === 1) return `${this.n}`;
    const sign = this.n < 0 ? "−" : "";
    const an = Math.abs(this.n);
    return `${sign}<sup>${an}</sup>&frasl;<sub>${this.d}</sub>`;
  }
  valueOf() { return this.n / this.d; }
}
