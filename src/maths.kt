/**
 * Created by say on 06.01.16.
 */
package maths
fun sin(v:Double): Double = Math.sin(v)
fun cos(v:Double): Double = Math.cos(v)
fun tan(v:Double): Double = Math.tan(v)
fun asin(v:Double): Double = Math.asin(v)
fun acos(v:Double): Double = Math.acos(v)
fun atan(v:Double): Double = Math.atan(v)
fun toRadians(v:Double): Double = Math.toRadians(v)
fun toDegrees(v:Double): Double = Math.toDegrees(v)
fun exp(v:Double): Double = Math.exp(v)
fun log(v:Double): Double = Math.log(v)
fun log10(v:Double): Double = Math.log10(v)
fun sqrt(v:Double): Double = Math.sqrt(v)
fun cbrt(v:Double): Double = Math.cbrt(v)
fun ceil(v:Double): Double = Math.ceil(v)
fun floor(v:Double): Double = Math.floor(v)
fun rint(v:Double): Double = Math.rint(v)
fun atan2(v:Double,x: Double): Double = Math.atan2(v, x)
fun pow(v:Double,exp: Double): Double = Math.pow(v, exp)
fun pow(v:Double,exp: Int): Double = Math.pow(v, exp.toDouble())
fun round(v:Double): Long = Math.round(v)
fun abs(v:Double): Double = Math.abs(v)
fun abs(v:Int):Int = Math.abs(v)
fun ulp(v:Double): Double = Math.ulp(v)
fun signum(v:Double): Double = Math.signum(v)
fun sinh(v:Double): Double = Math.sinh(v)
fun cosh(v:Double): Double = Math.cosh(v)
fun tanh(v:Double): Double = Math.tanh(v)
fun expm1(v:Double): Double = Math.expm1(v)
fun log1p(v:Double): Double = Math.log1p(v)
fun copySign(v:Double,sign: Double): Double = Math.copySign(v, sign)
fun exponent(v:Double): Int = Math.getExponent(v)
fun next(v:Double,direction: Double): Double = Math.nextAfter(v, direction)
fun nextUp(v:Double): Double = Math.nextUp(v)
fun scalb(v:Double,scaleFactor: Int): Double = Math.scalb(v, scaleFactor)
fun exp(c:Complex):Complex = c.exp()
fun pow(c:Complex,d:Complex) = c.pow(d)
fun pow(c:Complex,d:Double) = c.pow(d)
fun ln(c:Complex) = c.ln()
fun random():Double = Math.random()
operator fun Double.div(c:Complex): Complex = Complex(this,0.0)/c
operator fun Int.div(c:Complex): Complex = Complex(this.toDouble(),0.0)/c
operator fun Double.times(c:Complex): Complex = c*this
fun Double.toComplex() = Complex(this,0.0)
fun max(x:Double,y:Double):Double = if(x<y) y else x
fun max(x:Int,y:Int):Int = if(x<y) y else x
fun abs(x:Complex):Double = x.module()
fun real(x:Complex):Double = x.re
fun sqr(x:Double):Double = x*x
fun sqr(x:Int):Int = x*x
fun fact(n:Int):Int = if(n<=0) 1 else n*fact(n-1)

abstract class function() {
     abstract fun func(x: Double): Double
}


fun normmatrix(xx:Array<Array<Double>>):Double
{
    var norm = 0.0
    for(i in 0..xx.size-1)
        for(j in 0..xx[i].size-1)
            norm+=abs(xx[i][j])
    return norm/(xx.size)
}
fun normalizematrix(xx:Array<Array<Double>>):Array<Array<Double>>
{
    val nrm = normmatrix(xx)
    val c :Array<Array<Double>> = Array(xx.size,{i->Array(xx[i].size,{j->xx[i][j]/nrm})})
    return c
}

fun getvector(xx:Array<Array<Double>>,j:Int):Array<Double>
{
    val x:Array<Double> = Array(xx.size,{i->xx[i][j]})

    return x
}

fun printmatrix(x:Array<Array<Double>>)
{
    for(i in 0..x.size-1) {
        for (j in 0..x[i].size - 1)
            print("" + x[i][j] + "  ")
        println()
    }
    println()

}

fun lsgaus2(a:Array<Array<Double>>,b:Array<Array<Double>>):Array<Array<Double>>
{
    val n:Int  = a.size
    var f:Double
    var aa:Array<Array<Double>>
    var bb:Array<Array<Double>>
    var fl:Boolean
    var numfl:Double
    var maxa:Double
    val precise = 1e-24
    val x:Array<Array<Double>> = Array(n,{Array(n,{0.0})})
    maxa=abs(a[0][0])
    for(i in 0..n-1)
    for(j in 0..n-1)
        if(maxa<abs(a[i][j])) maxa=abs(a[i][j])

    fl=true
    numfl=1.0
   // do
   // {
        aa = Array(n,{i->Array(n,{j->a[i][j]/maxa})})
        bb = Array(n,{i->Array(n,{j->b[i][j]/maxa})})
        if(!fl) {
            for(i in 0..n-1)
            aa[i][i]=aa[i][i]+numfl*1e-10
        }
        fl=true


        for(i in 0..n-2)
        {

            if(abs(aa[i][i])<precise)
            {

                for(k in 0..n-1)
                if(k!=i)
                {
                    for(j in i+1..n-1)
                    aa[i][j]=aa[i][j]+aa[k][j]
                    for(j in 0..n-1)
                    bb[i][j]=bb[i][j]+bb[k][j]
                    if(abs(aa[i][i])>precise) break
                }
            }


            for(j in i+1..n-1)
            {
                f=aa[j][i]
                for(k in i..n-1)
                aa[j][k]=aa[j][k]-((aa[i][k]*f)/aa[i][i])


                for(k in 0..n-1)
                bb[j][k]=bb[j][k]-((bb[i][k]*f)/aa[i][i])

            }
        }

        for(k in 0..n-1) {
        x[n-1][k]=bb[n-1][k]/aa[n-1][n-1]
        for(i in n-2 downTo 0)
        {
            x[i][k]=bb[i][k];
            for(j in n-1 downTo i+1)
            x[i][k]=x[i][k]-aa[i][j]*x[j][k]
            x[i][k]=x[i][k]/aa[i][i]
        }
    }

      val m = mulmatrix(a,x)
        printmatrix(m)

        numfl=numfl*10
    //}
    //while(numfl<1e+10)
return x
}


fun inversematrix(a:Array<Array<Double>>):Array<Array<Double>>
{
    val n:Int = a.size
    val c:Array<Array<Double>> = Array(n,{Array(n,{0.0})})
    val v:Array<Array<Double>> = Array(n,{Array(n,{0.0})})
    var d:Double
    var b:Array<Array<Double>>
    for(j in 0..n-1)
    v[j][j]=1.0
    b= lsgaus2(a,v);


    for(i in 0..n-1)
    for(j in 0..n-1)
    {
        for(k in 0..n-1)
        c[i][j]=c[i][j]+b[i][k]*a[k][j]
    }

    d=0.0
    for(i in 0..n-1)
    for(j in 0..n-1)
    {
        if(i==j) d=d+abs(c[i][j]-1) else d=d+abs(c[i][j])
    }
    return b

}


fun mulmatrix(a:Array<Array<Double>>,b:Array<Array<Double>>):Array<Array<Double>>
{
    val n:Int = a.size
    val m:Int = a[0].size
    val k:Int = b[0].size
    val c:Array<Array<Double>> = Array(n,{Array(k,{0.0})})
    for(i in 0..n-1)
    for(j in 0..k-1)
    {
        for(l in 0..m-1)
        c[i][j]=c[i][j]+a[i][l]*b[l][j]
    }
    return c
}

fun mulvector(a:Array<Array<Double>>,b:Array<Double>):Array<Double>
{
    val n = a.size
    val m = b.size
    val c:Array<Double> = Array(n,{0.0})
    for(i in 0..n-1)
        for(j in 0..m-1)
           c[i]=c[i]+a[i][j]*b[j]
    return c
}

fun transp(a:Array<Array<Double>>):Array<Array<Double>>
{
    val  n:Int = a.size
    val  m:Int = a[0].size
    val b:Array<Array<Double>> = Array(m,{Array(n,{0.0})})


    for(i in 0..n-1)
    for(j in 0..m-1)
    b[j][i]=a[i][j]
    return b
}

fun integrate1(c:function,a:Double,b:Double,fa:Double,fb:Double,e:Double):Double
{

    //if((b-a)<e) return (c.func(a)+c.func(b))*(b-a)*0.5
    val val1=0.0
    val x = (b+a)/2.0
    val x1 = (x+a)/2.0
    val x2 = (x+b) /2.0

    val fx =  c.func(x)
    val fx1 =  c.func(x1)
    val fx2 =  c.func(x2)

    val i1 = (fa+fb)*0.5*(b-a)
    val dx2 = x-a
    val i2 = (fa+fb+fx*2)*dx2*0.5
    val dx3 = x2-x
    val i3 = (fa+2*fx1+2*fx+2*fx2+fb)*dx3*0.5
    val er = abs(i3-i1)+abs(i2-i1)+abs(i3-i2)
    if(er/(i3)<e) return i3
    else return integrate1(c,a,x1,fa,fx1,e)+integrate1(c,x1,x,fx1,fx,e)+
            integrate1(c,x,x2,fx,fx2,e)+integrate1(c,x2,b,fx2,fb,e)


}
fun integrate(c:function,a:Double,b:Double,e:Double):Double = integrate1(c,a,b,c.func(a),c.func(b),e)






class Complex(var re: Double, var im: Double) {


    // Constructors
    constructor(re: Double) : this(re, 0.0)
    constructor(re: Float) : this(re.toDouble(), 0.0)
    constructor(re: Int) : this(re.toDouble(), 0.0)
    constructor(re: Long) : this(re.toDouble(), 0.0)
    constructor(re: Short) : this(re.toDouble(), 0.0)
    // Unary operators
    fun exp(x:Complex):Complex = Complex(re, -im) // conjugate
    fun module():Double = sqrt(pow(re, 2) + pow(im, 2))
    // Comparison
    operator fun compareTo(that: Complex) = this.module().compareTo(that.module())

    // Arithmetic operations
    fun pow(c:Complex):Complex = (c*this.ln()).exp()
    fun pow(c:Double):Complex = (this.ln()*c).exp()

    fun exp():Complex  {
        val r: Double = exp(this.re)
        return Complex(r*cos(this.im),r*sin(this.im))
    }

    fun arg(): Double {
        return when  {
            this.re==0.0 && this.im>0.0 -> Math.PI/2
            this.re==0.0 && this.im<0.0 -> -Math.PI/2
            this.re>0.0  -> Math.atan(this.im/this.re)
            this.re<0.0 && this.im>=0.0 -> Math.atan(this.im/this.re)+Math.PI
            this.re<0.0 && this.im<0.0 -> Math.atan(this.im/this.re)-Math.PI
            else ->  0.0
        }
    }
    fun ln():Complex  {

        return Complex(log(this.module()),arg())
    }

    operator fun plus(c: Complex):Complex = Complex(re + c.re, im + c.im)
    operator fun minus(c: Complex):Complex = Complex(re-c.re,im-c.im)
    operator fun plus(c: Double):Complex = Complex(re + c, im)
    operator fun minus(c: Double):Complex = Complex(re - c, im)
    operator fun times(c: Complex):Complex =Complex(re * c.re - im * c.im, im * c.re + re * c.im)
    operator fun times(d: Double):Complex =Complex(re * d, im * d )
    operator fun div(c: Complex):Complex {
        val d = pow(c.re, 2) + pow(c.im, 2)
        return Complex((re * c.re + im * c.im) / d, (im * c.re - re * c.im) / d)
    }
    operator fun div(c: Double):Complex  = Complex(re / c, im / c)




    // String representation
    override fun toString() = asString()


    fun asString():String = this.re.toString() + (if (this.im < 0) "-" + -this.im else "+" + this.im) + "*i"



    // Factory methods
    fun apply(re: Double) = Complex(re)
    fun ln(c:Complex):Complex = c.ln()

    // Implicit conversions
    fun fromDouble(d: Double):Complex = Complex(d)
    fun fromFloat(f: Float):Complex =  Complex(f)
    fun fromLong(l: Long):Complex =  Complex(l)
    fun fromInt(i: Int):Complex =  Complex(i)
    fun fromShort(s: Short):Complex =  Complex(s)


}


fun Array<Double>.geti(x:Double):Int
{
    var i1:Int = 0
    var i2:Int = this.size-1

    while((i2-i1)>1)
    {
        val i = (i1+i2)/2
        if(x<this[i]) i2 = i else i1 = i
    }
    return i1

}

fun Array<Double>.reaxis(xx:Array<Double>,xn:Array<Double>):Array<Double>
{
    return Array<Double>(xn.size, {i->this.getx(xn[i],xx)})
}

fun Array<Double>.reaxis(xx:Array<Double>,a:Double,b:Double,n:Int):Array<Double>
{
    val h=(b-a)/n
    return Array<Double>(n, {i->this.getx(a+i*h,xx)})
}

fun Array<Double>.reaxis(xx:Array<Double>,a:Double,b:Double,h:Double):Array<Double>
{
    val n=((b-a)/h).toInt()
    return Array<Double>(n, {i->this.getx(a+i*h,xx)})
}

fun Array<Double>.smootmean(al:Int):Array<Double>
{
    val mas = Array<Double>(this.size,{0.0})
    for(i in al..this.size-1-al) {
        for (j in -al..al) {
            var i1 = i+j
            mas[i]+=this[i1]

        }
        mas[i] = mas[i]/(2*al+1)
    }
    for(i in 0..al-1) {
        for (j in -al..al) {
            var i1 = j + i
            if(i1<0) i1 = 0
            mas[i] += this[i1]

        }
        mas[i] = mas[i]/(2*al+1)
    }
    for(i in this.size-al..this.size-1) {
        for (j in -al..al) {
            var i1 = j + i
            if (i1 > this.size - 1) i1 = this.size - 1
            mas[i] += this[i1]

        }
        mas[i] = mas[i]/(2*al+1)
    }

    return mas
}

fun Array<Double>.getx(x:Double,xx:Array<Double>):Double
{
    val i = xx.geti(x)
    val a = (this[i+1]-this[i])/(xx[i+1]-xx[i])
    val b = this[i+1]-a*xx[i+1]
    return if(x<xx[0]) this[0] else if(x>xx.last()) this.last() else a*x+b
}

fun Array<Array<Double>>.getxy(x:Double,y:Double,xx:Array<Double>,yy:Array<Double>):Double
{
    val i = yy.geti(y)
    val v1 = this[i].getx(x,xx)
    val v2 = this[i+1].getx(x,xx)
    val a = (v2-v1)/(yy[i+1]-yy[i])
    val b = v2-a*yy[i+1]

    return if(y<yy.first()) v1 else if(y>yy.last()) v2 else a*y+b

}

fun Array<Array<Double>>.getxy(x:Double,y:Double,ax:Double,hx:Double,ay:Double,hy:Double):Double
{
    var i = ((y-ay)/hy).toInt()
    if(i<0) i = 0
    if(i>=this.size-1) i=this.size-2
    val v1 = this[i].getx(x,ax,hx)
    val v2 = this[i+1].getx(x,ax,hx)
    val a = (v2-v1)/(hy)
    val b = v1-a*(ay+hy*i)

    return if(y<=ay) v1 else if(y>=((this.size)*hy+ay)) v2 else a*y+b

}



fun Array<Double>.getx(x:Double,a:Double,h:Double):Double
{
    var i = ((x-a)/h).toInt()


    if(i<0) i = 0
    if(i>=this.size-1) i = this.size-2

    val aa = (this[i+1]-this[i])/h
    val b = this[i]-aa*(i*h+a)

    return (aa*x+b)
}

