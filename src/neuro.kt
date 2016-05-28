/**
 * Created by say on 08.04.16.
 */
package neuro
import maths.*

open class activation
{
    open fun act(g:Double):Double = 0.0
    open fun difact(g:Double):Double = 0.0
    open fun invact(g:Double):Double = 0.0

}

class logisticact:activation()
{
    override fun act(g:Double):Double = 1/(exp(-g)+1)
    override fun difact(g:Double):Double = sqr(act(g))*exp(-g)
    override fun invact(g:Double):Double =-log(1/g-1) //if(g<=0.001) -log(1/0.001-1) else if(g>=0.999) -log(1/0.999-1) else -log(1/g-1)
}

class neuron(cn:Int,caf:activation) {
    val af:activation = caf
    var n:Int = cn
    var w:Array<Double> = Array<Double>(n,{(random()-0.5)*0.01})
    var y = 0.0
    var s = 0.0
    fun ask(x:Array<Double>):Double {
        s=0.0
        for(i in 0..n-1)
            s = s+w[i]*x[i]
        y = af.act(s)
        return y
    }
}

class trainset(ntr:Int,cnin:Int,cnout:Int)
{
    var ntrain:Int = ntr
    var nin:Int = cnin
    var nout:Int = cnout
    var set:Array<Array<Array<Double>>> = Array(ntr,{arrayOf(Array(nin,{0.0}),Array(nout,{0.0}))})
    init {


    }


}

class inversetrain(nen:Int,nout:Int,nexmpl:Int,caf:activation)
{
    var w:Array<Array<Double>>  = Array(nout,{Array(nen,{(random()-0.5)*0.001})})
    val x:Array<Array<Double>> =  Array(nen,{Array(nexmpl,{0.0})})
    val af = caf
    var s1:Array<Double> = Array(nen,{0.0})
    var s2:Array<Double> = Array(nout,{0.0})
    val nexampl = nexmpl
    val nx = nen
    fun asky(x:Array<Double>):Array<Double>
    {
        val y = fvector(mulvector(w,x))
        return y

    }

    fun askxyx(x:Array<Double>):Array<Double>
    {
        s2 = mulvector(w,x)
        val y = fvector(s2)
        s1 = mulvector(transp(w),y)
        val xx = fvector(s1)
        return xx

    }


    fun ask(x:Array<Double>):Array<Double>
    {
        val d1 = mulvector(w,x)
        val y = fvector(d1)
        val d2 = mulvector(transp(w),y)
        val xx = fvector(d2)
        return xx

    }


    fun finvmatrix(x:Array<Array<Double>>):Array<Array<Double>>
    {
        val c:Array<Array<Double>> = Array(x.size,{i->Array(x[0].size,{j->af.invact(x[i][j])})})
        return c
    }

    fun fmatrix(x:Array<Array<Double>>):Array<Array<Double>>
    {
        val c:Array<Array<Double>> = Array(x.size,{i->Array(x[0].size,{j->af.act(x[i][j])})})
        return c
    }

    fun fvector(x:Array<Double>):Array<Double>
    {
        val c:Array<Double> = Array(x.size,{i->af.act(x[i])})
        return c
    }

    fun printmatr(x:Array<Array<Double>>)
    {
        for(i in 0..x.size-1) {
            for (j in 0..x[i].size - 1)
                print("" + x[i][j] + "  ")
        println()
        }
        println()

    }

    fun dataseterror():Double
    {
        var v=0.0
        for(i in 0..nexampl-1)
        {

            val xx=ask(getvector(x,i))
            for(j in 0..nx-1)
                v=v+abs(xx[j]-x[j][i])
        }

        return v

    }

    fun gradientbatch(al:Double,btch:Int)
    {
        for (i in 0..btch) {
            gradientlearn(al, (random() * nexampl).toInt())
            if((i%1000)==0) println(dataseterror())

        }
    }

    fun gradientlearn(al:Double,num:Int)
    {
        val dw:Array<Array<Double>> = Array(w.size,{i->Array(w[i].size,{0.0})})
        val xr = getvector(x,num)
        val xv = askxyx(xr)

        for(l in 0..w.size-1)
            for(p in 0..w[l].size-1) {
                for (k in 0..xr.size - 1)
                    dw[l][p] =dw[l][p] -(xr[k] - xv[k]) * af.difact(s1[k]) * af.difact(s2[l]) * w[l][k] * xr[p] - (xr[p] - xv[p]) * af.difact(s1[p]) * af.act(s2[l])

                w[l][p] =w[l][p]-al*dw[l][p]
            }
    }

    fun learnstep()
    {
        /*val xt = transp(x)
        val xx = inversematrix(mulmatrix(x,xt))
        val x1 = mulmatrix(xt,xx)
        val fix = finvmatrix(x)
        val wfx = mulmatrix(w,fix)
        val wi = inversematrix(mulmatrix(w,transp(w)))
        println(wi.size)
        println(wi[0].size)
        val ww1 = finvmatrix(mulmatrix(wi,wfx))
        printmatr(wi)


        w = mulmatrix(ww1,x1)

        val y = fmatrix(mulmatrix(w,x))
        val yt = transp(y)
        val p1= mulmatrix(fix,yt)
        val wt = mulmatrix(p1,inversematrix(mulmatrix(y,yt)))
        w = transp(wt)*/

        val fiy = mulmatrix(w,x)
        val xt = transp(x)
        val mxxt= mulmatrix(x,xt)
        val mxxti =inversematrix(mxxt)
        val xtx= mulmatrix(xt,mxxti)
        w = (mulmatrix(fiy,xtx))


        val y = fmatrix(mulmatrix(w,x))
        val yt = transp(y)
        val myyt = mulmatrix(y,yt)

        val myyti = inversematrix(myyt)


        val yty = mulmatrix(yt,myyti)
        val fix = finvmatrix(x)
        val wt = mulmatrix(fix,yty)
        val w = (transp(wt))
        //for(i in 0..ww.size-1)
        //    for(j in 0..ww[i].size-1)
        //        w[i][j] = (w[i][j]+ww[i][j])/2.0

        printmatr(w)

    }


}

open class neuronet(cnin:Int)
{
    lateinit var af:Array<activation>

    constructor(cnin:Int,n:Array<Int>,caf:Array<activation>): this(cnin)
    {
    setlayers(n,caf)
    af = caf
    }


    val nin:Int = cnin
    var nout:Int = 0
    lateinit var n:Array<Int>
    lateinit var net:Array<Array<neuron>>
    lateinit var y:Array<Array<Double>>

    fun setlayers(cn:Array<Int>,caf:Array<activation>)
    {
        n = cn
        net = Array<Array<neuron>>(n.size,{arrayOf()})
        net[0] = Array<neuron>(n[0],{neuron(nin,caf[0])})
        val l = n.size-1
        net[l] = Array<neuron>(n[l],{neuron(n[l-1],caf[l])})
        for(i in 1..l-1)
        {
            net[i] = Array<neuron>(n[i],{neuron(n[i-1],caf[i])})

        }
        nout = n[l]
        y = Array(n.size,{i->Array(n[i],{0.0})})

    }

    fun ask(x:Array<Double>):Array<Double>
    {

        if(x.size>nin) return arrayOf()
        for(i in 0..n[0]-1)
            y[0][i] = net[0][i].ask(x)

        for(j in 1..n.size-1)
            for(i in 0..n[j]-1)
                y[j][i] = net[j][i].ask(y[j-1])


        return y[n.size-1]
    }

    open fun learn()
    {

    }
    open fun learn(x:Array<Double>,y:Array<Double>)
    {}

}

open class backppg(cnin:Int,n:Array<Int>,caf:Array<activation>):neuronet(cnin,n,caf)
{
    val sg:Array<Array<Double>> = Array(n.size,{i->Array(n[i],{0.0})})
    var al = 0.01
    override fun learn(x:Array<Double>,y:Array<Double>)
    {

        var il=n.size-1
        var d = ask(x)
        for(i in 0..nout-1)
            sg[il][i] = (d[i]-y[i])*net[il][i].af.difact(net[il][i].s)


        il = n.size-2
        while(il>=0) {
            for(i in 0..n[il]-1)
            {
                sg[il][i] = 0.0
                for(k in 0..n[il+1]-1) {


                    sg[il][i] = sg[il][i] + sg[il + 1][k] * net[il + 1][k].w[i]
                }
                sg[il][i] = sg[il][i] * net[il][i].af.difact(net[il][i].s)
            }
            il = il - 1
        }
        il = 0

        for(i in 0..n[il]-1) {
            var al1 = 0.0
            if (random() > 0.999999)
                al1 = 0.01
            for(k in 0..net[il][i].n-1)
                net[il][i].w[k] = net[il][i].w[k] - al * sg[il][i] * x[k] + al1 * (random() - 0.5)
        }

        for(il in 1..n.size-1) {
            for(i in 0..n[il]-1) {
                var al1 = 0.0
                if (random() > 0.999999)
                    al1 = 0.1
                for(k in 0..net[il][i].n-1)
                    net[il][i].w[k] = net[il][i].w[k] - al * sg[il][i] * net[il - 1][k].y + al1 * (random() - 0.5)
            }
        }


    }


}

class prelearngradient(cnin:Int,n:Array<Int>,caf:Array<activation>,l:trainset):backppg(cnin,n,caf)
{
    lateinit var encoders:Array<inversetrain>
    init {

        encoders = arrayOf(inversetrain(nin+nout,n[0],l.ntrain,af[0]))+Array(n.size-1,{i->inversetrain(n[i],n[i+1],l.ntrain,af[i])})

        for(i in 0..l.ntrain-1) {
            for (j in 0..nin - 1)
                encoders[0].x[j][i] = l.set[i][0][j]
            for (j in nin..nin + nout - 1)
                encoders[0].x[j][i] = l.set[i][1][j-nin]
        }

        encoders[0].gradientbatch(0.1,10000)

        var ans:Array<Array<Double>> = Array(l.ntrain,{Array(0,{0.0})})
        for(i in 0..l.ntrain-1)
          ans[i]=encoders[0].asky(getvector(encoders[0].x,i))

        for(j in 1..encoders.size-1)
        {

            println(encoders[j].nx)
            println(ans.size)
            println(ans[0].size)
            for (i in 0..encoders[j].nx - 1)
                for (k in 0..encoders[j].nexampl - 1)
                encoders[j].x[i][k]= ans[k][i]
               encoders[j].gradientbatch(0.1,10000)
               for(i in 0..l.ntrain-1)
               ans[i]=encoders[j].asky(getvector(encoders[j].x,i))
        }

        for(i in 0..n.size-2)
        {
            for (j in 0..n[i]-1)
                for (k in 0..net[i][j].n-1)
            net[i][j].w[k] = encoders[i].w[j][k]

        }


    }

}
