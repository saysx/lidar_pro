/**
 * Created by say on 07.01.16.
 */
package mie

import maths.*
import java.io.*
import java.sql.Connection
import java.sql.DriverManager
import java.sql.Statement
import java.util.*
data class optpar(val n:Int)
{
    var sct:Double = 0.0
    var ext:Double = 0.0
    var spi:Double = 0.0
    val ind: Array<Double> = Array(n,{i->0.0})
}




fun mie(refrp:Complex,refrm:Complex,nang:Int,ang:Array<Double>,r:Double,lam:Double):optpar
{
    val opt : optpar = optpar(nang)
    val x:Double = 2.0*Math.PI*r/lam
    val nmax:Int = (x+4.0*pow(x,1.0/3.0)+2.0).toInt()
    val m:Complex = refrp/refrm
    val mx:Complex = x*m
    val nd:Int = max(nmax,abs(mx).toInt())+25
    val D:Array<Complex> = Array(nd+1,{i->(i+1).toDouble()/mx})
    D[nd] = Complex(0.0,0.0)
    for(i in nd downTo 1) {
        D[i - 1] = D[i - 1] - 1.0 / (D[i] + D[i - 1])
    }
    val psi = Array<Double>(nd+1,{0.0})
    val ksi = Array<Double>(nd+1,{0.0})
    val eps : Array<Complex>
    psi[0]=sin(x)
    psi[1]=sin(x)/x-cos(x)
    ksi[0]=cos(x)
    ksi[1]=cos(x)/x+sin(x)

    for(i in 2..nd)
    {
        psi[i]=((2.0*(i-1)+1.0)*psi[i-1]/x)-psi[i-2]
        ksi[i]=((2.0*(i-1)+1.0)*ksi[i-1]/x)-ksi[i-2]
    }

    eps = Array<Complex>(nd+1,{i->Complex(psi[i],-ksi[i])})
    val a = Array<Complex>(nd+1,{Complex(0.0,0.0)})
    val b = Array<Complex>(nd+1,{Complex(0.0,0.0)})
    for(i in 1..nd)
    {
        val tmp1=D[i]/m+i/x
        val tmp2=D[i]*m+i/x
        a[i] = (tmp1 * psi[i] - psi[i - 1]) / (tmp1 * eps[i] - eps[i - 1])
        b[i] = (tmp2 * psi[i] - psi[i - 1]) / (tmp2 * eps[i] - eps[i - 1])
    }

    val nang1 = 45
    val pq1 =Array<Double>(nd+1,{0.0})
    val tq1 =Array<Double>(nd+1,{0.0})
    val s11 = Array<Complex>(nang1+1,{Complex(0.0,0.0)})
    val s21 = Array<Complex>(nang1+1,{Complex(0.0,0.0)})

    for(j in 0..nang1-1)
    {
        pq1[0]=0.0
        pq1[1]=1.0
        tq1[1]=cos(Math.PI*j/nang1)
        s11[j] = Complex(0.0,0.0)
        s21[j] = Complex(0.0,0.0)
        for(i in 2..nd)
        {
            pq1[i]=((2.0*i-1.0)/(i-1.0))*tq1[1]*pq1[i-1]-i*pq1[i-2]/(i-1.0);
            tq1[i]=i*tq1[1]*pq1[i]-(i+1.0)*pq1[i-1];
        }

        for(i in 1..nd-1)
        {
            s11[j]=s11[j]+(2.0*i+1.0)*(a[i]*pq1[i]+b[i]*tq1[i])/((i+1.0)*i)
            s21[j]=s21[j]+(2.0*i+1.0)*(a[i]*tq1[i]+b[i]*pq1[i])/((i+1.0)*i)
        }
    }

    var indi:Double = 0.0
    for(i in 0..nang1-1)
      indi = indi + (pow(abs(s11[i]),2)+pow(abs(s21[i]),2)+pow(abs(s11[i+1]),2)+pow(abs(s21[i+1]),2))*0.5*Math.PI/nang1;

    val pq = Array<Double>(nd+1,{0.0})
    val tq = Array<Double>(nd+1,{0.0})
    val s1 = Array<Complex>(nang,{Complex(0.0,0.0)})
    val s2 = Array<Complex>(nang,{Complex(0.0,0.0)})

    for(j in 0..nang-1)
    {
        pq[0]=0.0;
        pq[1]=1.0;
        tq[1]=cos(ang[j]);
        s1[j] = Complex(0.0,0.0);
        s2[j] = Complex(0.0,0.0);
        for(i in 2..nd)
        {
            pq[i]=((2.0*i-1.0)/(i-1.0))*tq[1]*pq[i-1]-i*pq[i-2]/(i-1.0)
            tq[i]=i*tq[1]*pq[i]-(i+1.0)*pq[i-1]
        }
        for(i in 1..nd)
        {
            s1[j]=s1[j]+(2.0*i+1.0)*(a[i]*pq[i]+b[i]*tq[i])/((i+1.0)*i);
            s2[j]=s2[j]+(2.0*i+1.0)*(a[i]*tq[i]+b[i]*pq[i])/((i+1.0)*i);
        }
    }

    var qsct=0.0
    var qext=0.0
    for(i in 1..nd)
    {
        qext = qext+(2.0*i+1.0)*real(a[i]+b[i])
        qsct = qsct+(2.0*i+1.0)*(pow(abs(a[i]),2)+pow(abs(b[i]),2))
    }
    qext=qext*2.0/pow(x,2.0);
    qsct=qsct*2.0/pow(x,2.0);
    opt.sct=qsct;
    opt.ext=qext;
    opt.spi = pow(abs(s1[nang-1]),2)/pow(x,2);
    opt.spi=opt.spi/Math.PI;
    for(i in 0..nang-1)
        opt.ind[i] = (pow(abs(s1[i]),2)+pow(abs(s2[i]),2))/indi;

    return opt
}
class hitran
{


}

object Ph {
    val Kbolc:Double = 1.3806485279e-23
    fun bar_to_at(p:Double):Double = p*1.019716
    fun bar_to_Pa(p:Double):Double = p*1.0e+5
    fun at_to_bar(p:Double):Double = p*0.980665205
    fun mbar_to_bar(p:Double):Double = p*1e-3
    fun cm_1_to_mkm(lam:Double):Double = 10000/lam
    fun mkm_to_cm_1(wave:Double):Double = 10000/wave
    fun Pa_to_atm(p:Double):Double = 9.86923e-6*p


    //p - pressure [Pa],T - temperature [K], lam - wavelength [mkm]
    fun nair(p:Double,T:Double,lam:Double):Double =  1.0+p*1e-8*(77.6+0.584/pow(lam,2))/T
    // N-concentration in [cm-3] T - temperature[K] P- pressure in [Pa]
    fun NinP_air(N:Double,T:Double):Double = N*Kbolc*T*1e+6
    // P - pressure in [Pa] result N in [cm-3]
    fun PinN_air(P:Double,T:Double):Double = P*1e-6/(Kbolc*T)

}


open class GaseModel
{
    //[atm]
    open fun getp(x:Double,y:Double,z:Double):Double =1.0
    //[K]
    open fun gett(x:Double,y:Double,z:Double):Double = 1.0
    //[atm]
    open fun getg(num:Int,x:Double,y:Double,z:Double):Double = 1.0
    open fun getnumbyname(name:String):Int = 1
    open fun maxgases():Int = 1
    open fun getnamebynum(num:Int):String = ""
    //[km]
    open fun maxz():Double = 100.0
    open fun maxx():Double = 0.0
    open fun maxy():Double = 0.0
    open fun minz():Double = 0.0
    open fun minx():Double = 0.0
    open fun miny():Double = 0.0

}

class GaseModelUsa(val fn:String) : GaseModel()
{
     init {
          readUSAMet(fn)
          }




    val usagase:Array<String> = arrayOf("H2O","CO2","O3","N2O","CO","CH4","O2","NO","SO2","NO2","NH3","HNO3","OH","HF","HCl", "HBr", "HI", "ClO", "OCS", "H2CO", "HOCl", "N2", "HCN", "CH3Cl", "H2O2", "C2H2", "C2H6", "PH3", "COF2", "SF6", "H2S", "HCOOH", "HO2", "O", "ClONO2", "NO+", "HOBr", "C2H4")

    override fun maxgases():Int = ng-1
    override fun getnamebynum(num:Int):String = if (num>=0 && num<maxgases()) usagase[num] else ""


    override fun getp(x:Double,y:Double,z:Double):Double =  p.getx(z,a,hh)

    override fun gett(x:Double,y:Double,z:Double):Double = t.getx(z,a,hh)
    override fun getg(num:Int,x:Double,y:Double,z:Double):Double =  g[num].getx(z,a,hh)
    override fun getnumbyname(name:String):Int = usagase.indexOf(name)
    override fun maxz():Double = h[h.size-1]
    override fun minz():Double = 0.0

    public var nh:Int = 0
    private var ng:Int = 0
    lateinit private var h:Array<Double>
    lateinit private var g:Array<Array<Double>>
    lateinit private var p:Array<Double>
    lateinit private var t:Array<Double>
    private var a:Double = 0.0
    private var b:Double = 0.0
    private var hh:Double = 0.0
    private var n:Int = 0

    fun readUSAMet(fn:String) {

        val sc = Scanner(File(fn));

        sc.useLocale(Locale.US)

        nh = sc.nextInt()
        ng = sc.nextInt()
        h = Array<Double>(nh, { 0.0 })
        t = Array<Double>(nh, { 0.0 })
        p = Array<Double>(nh, { 0.0 })
        val proc = Array<Double>(nh, { 0.0 })
        g = Array(ng, { Array<Double>(nh, { 0.0 }) })

        for (i in 0..30)
            sc.next()


        for (i in 0..nh - 1) {
            h[i] = sc.nextDouble()

            p[i] = sc.nextDouble()*100
            t[i] = sc.nextDouble()


            for (j in 0..ng - 1) {
                g[j][i] = sc.nextDouble()
                proc[i] = proc[i]+g[j][i]
            }

        }

        for (i in 0..nh - 1)
            for (j in 0..ng - 1)
                g[j][i] = g[j][i]*p[i]/proc[i]

        a = h[0]
        b = h[nh - 1]
        hh = 0.1
        n = ((b - a) / hh).toInt()
        p = p.reaxis(h, a, b, hh)
        t = t.reaxis(h, a, b, hh)
        for (j in 0..ng - 1)
            g[j] = g[j].reaxis(h, a, b, hh)


        println("sx = " + t.size)
    }





}

open class AerosolModel()
{
    open fun getext(lam:Double,x:Double,y:Double,z:Double):Double = 0.0 //lam - mkm,x,y,z - km, [km-1]
    open fun getsct(lam:Double,x:Double,y:Double,z:Double):Double = 0.0
    open fun getabs(lam:Double,x:Double,y:Double,z:Double):Double = 0.0
    open fun getsctang(lam:Double,x:Double,y:Double,z:Double,al:Double):Double = 0.0
    open fun getind(lam:Double,x:Double,y:Double,z:Double,al:Double):Double = 0.0
}

class AerosolModelVar:AerosolModel() {

    init {
        readdb()

    }

    val smalldisp: Array<Double> = arrayOf(1.0, 0.0, 41.0, 26.0, 5.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 10.0, 0.0, 0.0)
    val accumdisp: Array<Double> = arrayOf(5.0, 0.0, 34.0, 19.0, 0.0, 0.0, 0.0, 3.0, 2.0, 14.0, 0.0, 13.0, 10.0, 0.0, 0.0)
    val roughdisp: Array<Double> = arrayOf(5.0, 0.0, 15.0, 0.0, 1.0, 0.0, 0.0, 10.0, 23.0, 0.0, 0.0, 35.0, 0.0, 0.0, 0.0)
    val density: Array<Double> = arrayOf(1.0, 0.9, 1.77, 2.55, 2.0, 0.0, 0.0, 4.0, 5.3, 2.17, 0.0, 2.65, 1.05, 1.25, 1.25)
    val nc: Array<Double> = arrayOf(20000.0, 200.0, 2.5)
    lateinit var lam: Array<Array<Double>>
    lateinit var refre: Array<Array<Double>>
    lateinit var refim: Array<Array<Double>>

    fun makeconc(h: Array<Double>, conc: Array<Double>)
    {
    
    }

    fun getextsct(lam:Double,h:Array<Double>,ext:Array<Double>,sct:Array<Double>)
    {

    }

    fun readdb()
    {
        Class.forName("org.sqlite.JDBC")

        var conn: Connection
        var state: Statement


        conn = DriverManager.getConnection("jdbc:sqlite:common.db")
        state = conn.createStatement()
        conn.setAutoCommit(false)

        var res = state.executeQuery("select count(DISTINCT num) as cnt from REFR")
        res.next()
        val numf = res.getInt("cnt")
        println(numf)
        lam = Array(numf,{Array(0,{0.0})})
        refre = Array(numf,{Array(0,{0.0})})
        refim = Array(numf,{Array(0,{0.0})})

        for(i in 0..numf-1)
        {
            res = state.executeQuery("select count(lam) as cnt from REFR where num="+i+"")
            res.next()
            val nlam = res.getInt("cnt")
            lam[i] = Array(nlam,{0.0})
            refre[i] = Array(nlam,{0.0})
            refim[i] = Array(nlam,{0.0})

            var res = state.executeQuery("select lam,refrre,refrim from REFR where num="+i+" order by lam")
            var j=0
            while(res.next()||j<nlam)
            {
                lam[i][j] = res.getDouble("lam")
                refre[i][j] = res.getDouble("refrre")
                refim[i][j] = res.getDouble("refrim")
                println(" "+lam[i][j]+" "+refre[i][j]+" "+refim[i][j])
                j++


            }


        }

    }

    fun crefr(m:Array<Double>,proc:Array<Double>):Double
    {
        var c = 0.0
        var s = 0.0
        for(i in 0..m.size-1) {
            s = s + proc[i] / density[i]
            c = c + m[i] * proc[i] / density[i]
        }
        return c / s
    }

    override fun getext(lam1:Double,x:Double,y:Double,z:Double):Double = 0.0 //lam - mkm,x,y,z - km, [km-1]
    override fun getsct(lam:Double,x:Double,y:Double,z:Double):Double = 0.0
    override fun getabs(lam:Double,x:Double,y:Double,z:Double):Double = 0.0


}


object AerosolModelIAO:AerosolModel()
{

    override fun getext(lam1:Double,x:Double,y:Double,z:Double):Double = ext.getxy(z,lam1,h,lam) //lam - mkm,x,y,z - km, [km-1]
    override fun getsct(lam1:Double,x:Double,y:Double,z:Double):Double = 0.0
    override fun getabs(lam1:Double,x:Double,y:Double,z:Double):Double = 0.0
    override fun getsctang(lam1:Double,x:Double,y:Double,z:Double,al:Double):Double = sctpi.getxy(z,lam1,h,lam)

    private lateinit var lam:Array<Double>
    private lateinit var h:Array<Double>
    private lateinit var ext:Array<Array<Double>>
    private lateinit var sctpi:Array<Array<Double>>
    init {
        Class.forName("org.sqlite.JDBC")

        var conn: Connection
        var state: Statement


        conn = DriverManager.getConnection("jdbc:sqlite:common.db")
        state = conn.createStatement()
        conn.setAutoCommit(false)

        var res = state.executeQuery("select count(DISTINCT lam) as cnt from AEREXT")
        res.next()
        val nl = res.getInt("cnt")
        res = state.executeQuery("select count(DISTINCT height) as cnt from AEREXT")
        res.next()
        val nh = res.getInt("cnt")

        res = state.executeQuery("select value from AEREXT ORDER BY lam, height")

        ext = Array(nl,{Array<Double>(nh,{0.0})})
        for(i in 0..nl-1)
            for(j in 0..nh-1)
            {
                res.next()
                ext[i][j] = res.getDouble("value")

            }


        res = state.executeQuery("select value from AERSCTPI ORDER BY lam, height")

        sctpi = Array(nl,{Array<Double>(nh,{0.0})})
        for(i in 0..nl-1)
            for(j in 0..nh-1)
            {
                res.next()
                sctpi[i][j] = res.getDouble("value")

            }


        res = state.executeQuery("select distinct lam from AEREXT ORDER BY lam")
        lam = Array<Double>(nl,{0.0})
        for(i in 0..nl-1)
        {
            res.next()
            lam[i] = res.getDouble("lam")
        }

        res = state.executeQuery("select distinct height from AEREXT ORDER BY height")
        h = Array<Double>(nh,{0.0})
        for(i in 0..nh-1)
        {
            res.next()
            h[i] = res.getDouble("height")
        }



        res = state.executeQuery("select value from AEREXT ORDER BY lam, height")

        ext = Array(nl,{Array<Double>(nh,{0.0})})
        for(i in 0..nl-1)
            for(j in 0..nh-1)
            {
                res.next()
                ext[i][j] = res.getDouble("value")

            }

        res = state.executeQuery("select distinct lam from AEREXT ORDER BY lam")
        lam = Array<Double>(nl,{0.0})
        for(i in 0..nl-1)
        {
            res.next()
            lam[i] = res.getDouble("lam")
        }

        res = state.executeQuery("select distinct height from AEREXT ORDER BY height")
        h = Array<Double>(nh,{0.0})
        for(i in 0..nh-1)
        {
            res.next()
            h[i] = res.getDouble("height")
        }



        conn.commit()
        state.close()
        conn.close()


    }


}

object GaseModelDbUsa:GaseModel()
{



    val usagase:Array<String> = arrayOf("H2O","CO2","O3","N2O","CO","CH4","O2","NO","SO2","NO2","NH3","HNO3","OH","HF","HCl", "HBr", "HI", "ClO", "OCS", "H2CO", "HOCl", "N2", "HCN", "CH3Cl", "H2O2", "C2H2", "C2H6", "PH3", "COF2", "SF6", "H2S", "HCOOH", "HO2", "O", "ClONO2", "NO+", "HOBr", "C2H4")

    override fun maxgases():Int = ng-1
    override fun getnamebynum(num:Int):String = if (num>=0 && num<maxgases()) usagase[num] else ""


    override fun getp(x:Double,y:Double,z:Double):Double =  p.getx(z,a,hh)

    override fun gett(x:Double,y:Double,z:Double):Double = t.getx(z,a,hh)
    override fun getg(num:Int,x:Double,y:Double,z:Double):Double =  g[num].getx(z,a,hh)
    override fun getnumbyname(name:String):Int = usagase.indexOf(name)
    override fun maxz():Double = h[h.size-1]
    override fun minz():Double = 0.0


    private var ng:Int = 0
    lateinit private var h:Array<Double>
    lateinit private var g:Array<Array<Double>>
    lateinit private var p:Array<Double>
    private var t:Array<Double>
    private var a:Double = 0.0
    private var b:Double = 0.0
    private var hh:Double = 0.0
    private var n:Int = 0




    init {
        Class.forName("org.sqlite.JDBC")

        var conn: Connection
        var state: Statement


            conn = DriverManager.getConnection("jdbc:sqlite:common.db")
            state = conn.createStatement()
            conn.setAutoCommit(false)


                val res1 = state.executeQuery("SELECT count(*) as cnt FROM sqlite_master WHERE type='table' AND name='USA_MODEL'")
                res1.next()
                if (res1.getInt("cnt") > 0) {
                    println("Ok start db")
                } else {

                    conn.close()
                    val gm = CommonDbModel()
                    conn = DriverManager.getConnection("jdbc:sqlite:common.db")
                    state = conn.createStatement()
                    conn.setAutoCommit(false)
                    println("reinsert data to db")
                }


        var res = state.executeQuery("select count(*) as cnt from USA_MODEL where parid=-1")
        res.next()
        n = res.getInt("cnt")

        res = state.executeQuery("select count(*) as cnt from USA_PARS where parid>=0")
        res.next()
        ng = res.getInt("cnt")

        println(ng)

        t = Array<Double>(n,{0.0})
        p = Array<Double>(n,{0.0})
        g = Array(n,{Array<Double>(n,{0.0})})

        res = state.executeQuery("select value from USA_MODEL where parid=-1")
        var i =0
        while(res.next()) {
            t[i] = res.getDouble("value")
            i++
            if(i>=n) break
        }

        res = state.executeQuery("select value from USA_MODEL where parid=-2")
        i =0
        while(res.next()) {
            p[i] = res.getDouble("value")
            i++
            if(i>=n) break
        }
        for(j in 0..ng-1)
        {
            res = state.executeQuery("select value from USA_MODEL where parid="+j)
            i =0
            while(res.next()) {
                g[j][i] = res.getDouble("value")
                i++
                if(i>=n) break
            }
        }
        res = state.executeQuery("select max(height) from USA_MODEL where parid=0")
        res.next()
        b = res.getDouble(1)
        res = state.executeQuery("select min(height) from USA_MODEL where parid=0")
        res.next()
        a = res.getDouble(1)
        hh = (b-a)/n
        println(a)
        println(b)
        conn.commit()
        state.close()
        conn.close()




    }

}






class CommonDbModel()
{

fun addaer()
{

    Class.forName("org.sqlite.JDBC")
    val conn: Connection = DriverManager.getConnection("jdbc:sqlite:common.db")
    val state: Statement = conn.createStatement()
    conn.setAutoCommit(false)

    state.execute("DROP TABLE IF EXISTS AEREXT")
    state.execute("DROP TABLE IF EXISTS AERSCTPI")

    state.execute("CREATE TABLE AEREXT ("+
            "lam float, height float, value float,"+
            "primary key(lam,height))")
    state.execute("CREATE TABLE AERSCTPI ("+
            "lam float, height float, value float,"+
            "primary key(lam,height))")

    val f = File("dat/AEROZOL.DAT")
    val sc= Scanner(f)
    sc.useLocale(Locale.US)
    val nl = sc.nextInt()
    val nh = sc.nextInt()

    val lams = Array<Double>(nl,{0.0})
    for(i in 0..nl-1)
        lams[i] = sc.nextDouble()


    val sql1:String = "insert into  AEREXT values(?,?,?)"
    val preps1 = conn.prepareStatement(sql1)
    for(j in 0..nh-1) {
        val h = sc.nextDouble()
        for (i in 0..nl - 1) {
            val ext = sc.nextDouble()
            preps1.setDouble(1, lams[i])
            preps1.setDouble(2, h)
            preps1.setDouble(3, ext)
            preps1.addBatch()

        }
    }

        val sql2:String = "insert into  AERSCTPI values(?,?,?)"
        val preps2 = conn.prepareStatement(sql2)
        for(j in 0..nh-1) {
            val h = sc.nextDouble()
            for (i in 0..nl - 1)
            {
                val sctpi = sc.nextDouble()
                preps2.setDouble(1,lams[i])
                preps2.setDouble(2,h)
                preps2.setDouble(3,sctpi)
                preps2.addBatch()

            }
        }


    preps1.executeBatch()
    preps2.executeBatch()
    state.close()
    conn.commit()

    conn.close()


}

fun addstatsum() {
    Class.forName("org.sqlite.JDBC")
    val conn: Connection = DriverManager.getConnection("jdbc:sqlite:common.db")
    val state: Statement = conn.createStatement()
    conn.setAutoCommit(false)

    state.execute("DROP TABLE IF EXISTS STATSUM")


    state.execute("CREATE TABLE STATSUM ("+
            "num int, temp float, value float,"+
            "primary key(num,temp))")

    /*val f = File("dat/")
    val files = f.listFiles(object: FilenameFilter {
        override fun accept(dir: File, name: String): Boolean =
                if (name.indexOf("stat.dat") > 0)  true
                else false
    })*/


    val statf:Array<String> = arrayOf("ph3stat.dat", "c2h6stat.dat", "c2h2stat.dat", "h2o2stat.dat",
        "ch3clstat.dat","hcnstat.dat","n2stat.dat", "hoclstat.dat", "h2costat.dat", "ocsstat.dat",
        "clostat.dat", "histat.dat", "hbrstat.dat","hclstat.dat", "hfstat.dat","ohstat.dat",
        "hno3stat.dat", "nh3stat.dat", "no2stat.dat","so2stat.dat", "nostat.dat", "o2stat.dat",
        "ch4stat.dat", "costat.dat", "n2ostat.dat", "o3stat.dat", "co2stat.dat", "h2ostat.dat")
    val statf1:Array<String> = Array<String>(28,{""})
    for(i in 0..27)
        statf1[i] = statf[27-i]



    val strs = Array<String>(28,{""})
    val sql1:String = "insert into  STATSUM values(?,?,?)"
    val preps1 = conn.prepareStatement(sql1)


    for(i in 0..27)
    {
        val f = File("dat/"+statf1[i])
        val sc= Scanner(f)
        sc.useLocale(Locale.US)
        strs[i] = statf1[i].substringBefore("stat.dat").toUpperCase()
        println(strs[i])

        while(sc.hasNextDouble())
        {
            val temp = sc.nextDouble()
            val value = sc.nextDouble()


            preps1.setInt(1,i)
            preps1.setDouble(2,temp)
            preps1.setDouble(3,value)
            preps1.addBatch()
        }


    }
    preps1.executeBatch()
    state.close()
    conn.commit()

    conn.close()

}

fun addrefr()  {
    Class.forName("org.sqlite.JDBC")
    val conn: Connection = DriverManager.getConnection("jdbc:sqlite:common.db")
    val state: Statement = conn.createStatement()
    conn.setAutoCommit(false)

    state.execute("DROP TABLE IF EXISTS REFR")


    state.execute("CREATE TABLE REFR ("+
    "num int, lam float, refrre float, refrim float,"+
    "primary key(num,lam))")

    val NameF = arrayOf("water.dat","iceT273.dat", "amonsulf.dat","carbon.dat","soot.dat","meteordst.dat","volcdust.dat",
            "al2o3all.dat","hematiteT293.dat","seasult.dat","ariddst.dat","quarz.dat","organic.dat",
            "h2s04215.dat", "h2s04315.dat", "ethglyc.dat", "SiO2_kf.dat", "al2o3_kr.dat")


    val sql1:String = "insert into  REFR values(?,?,?,?)"
    val preps1 = conn.prepareStatement(sql1)

    for(i in 0..NameF.size-1)
    {
        val f = File("dat/"+NameF[i])
        val sc= Scanner(f)
        sc.useLocale(Locale.US)
        while(sc.hasNext())
        {
            val lam = sc.nextDouble()
            val refre = sc.nextDouble()
            val refim = sc.nextDouble()

            preps1.setInt(1,i)
            preps1.setDouble(2,lam)
            preps1.setDouble(3,refre)
            preps1.setDouble(4,refim)
            preps1.addBatch()


        }
    }
    preps1.executeBatch()
    state.close()
    conn.commit()

    conn.close()

}


fun    addxsc() {
    Class.forName("org.sqlite.JDBC")
    val conn: Connection = DriverManager.getConnection("jdbc:sqlite:common.db")
    val state: Statement = conn.createStatement()
    conn.setAutoCommit(false)

    state.execute("DROP TABLE IF EXISTS XSC")
    state.execute("DROP TABLE IF EXISTS XSCPAR")

    state.execute("CREATE TABLE IF NOT EXISTS XSCPAR (" +
    "id int primary key,"+
    "mol varchar(7),"+
    "temperature float,"+
    "startlam float, endlam float)"
    )

    state.execute("CREATE TABLE IF NOT EXISTS XSC (" +
            "id int, idval int,"+

            "value float, primary key(id,idval), foreign key (id) references XSCPAR(id))"
    )

    val f = File("dat/")
    val files = f.listFiles(object: FilenameFilter {
        override fun accept(dir: File, name: String): Boolean =
                if (name.indexOf(".xsc") > 0)  true   //если встречаем в имени файла ".rtf"
                else false
    })

    val sql1:String = "insert into  XSCPAR values(?,?,?,?,?)"
    val preps1 = conn.prepareStatement(sql1)
    val sql2:String = "insert into  XSC values(?,?,?)"
    val preps2 = conn.prepareStatement(sql2)


    for(i in 0..files.count()-1)
    {



        val sc= Scanner(files[i])
        sc.useLocale(Locale.US)
        val name = sc.findWithinHorizon(".{20}",20).trim()
        print(name+" ")
        val minwave = sc.findWithinHorizon(".{10}",10).trim().toDouble()
        print(minwave.toString()+" ")
        val maxwave = sc.findWithinHorizon(".{10}",10).trim().toDouble()
        print(maxwave.toString()+" ")
        val numpoint = sc.findWithinHorizon(".{7}",7).trim().toInt()
        print(numpoint.toString()+" ")

        val temp = sc.findWithinHorizon(".{7}",7).trim().toDouble()
        print(temp.toString()+" ")
        val press = sc.findWithinHorizon(".{6}",6).trim().toDouble()
        print(press.toString()+" ")
        val maxcs = sc.findWithinHorizon(".{10}",10).trim().toDouble()
        print(maxcs.toString()+" ")
        val instrres = sc.findWithinHorizon(".{5}",5).trim()
        print(instrres.toString()+" ")
        val cname = sc.findWithinHorizon(".{15}",15).trim()
        print("cname= "+cname+" ")
        val res1 = sc.findWithinHorizon(".{4}",4).trim()
        print("res1= "+res1+" ")
        val res2 = sc.findWithinHorizon(".{3}",3).trim()
        print("res2= "+ res2+" ")
        val ref = sc.findWithinHorizon(".{3}",3).trim()
        print("ref = "+ref.toString()+" ")
        println()
        preps1.setInt(1,i)
        preps1.setNString(2,name)
        preps1.setDouble(3,temp)
        preps1.setDouble(4,minwave)
        preps1.setDouble(5,maxwave)
        preps1.addBatch()

        for(j in 0..numpoint-1)
        {
            val val1 = sc.nextDouble()
            preps2.setInt(1,i)
            preps2.setInt(1,j)
            preps2.setDouble(1,val1)
            preps2.addBatch()
        }




    }

    preps1.executeBatch()
    preps2.executeBatch()
    conn.commit()
    conn.close()

}

    init {

     Class.forName("org.sqlite.JDBC")

     val conn: Connection = DriverManager.getConnection("jdbc:sqlite:common.db")
     val state: Statement = conn.createStatement()
     conn.setAutoCommit(false)




     state.execute("DROP TABLE IF EXISTS USA_MODEL")
     state.execute("DROP TABLE IF EXISTS USA_PARS")

     state.execute("CREATE TABLE IF NOT EXISTS USA_PARS (" +
             "parid Int primary key,"+
             "parname varchar(15),"+

             "parmes varchar(5)"+
             ")")


     state.execute("CREATE TABLE IF NOT EXISTS USA_MODEL (" +
             "height float," +
             "parid Int,"+
             "value float," +
             "primary key(parid,height),"+
             "foreign key (parid) references USA_PARS(parid)"+
                          ")")


     val gm = GaseModelUsa("dat/usa.met")
     val h:Double = 0.1
     val nh:Int = (gm.maxz()/h).toInt()


     state.execute("insert into  USA_PARS values(-2,'temperature','K')")
     state.execute("insert into  USA_PARS values(-1,'pressure','Pa')")
     for(i in 0..gm.maxgases())
         state.execute("insert into  USA_PARS values("+i.toString()+",'"+gm.getnamebynum(i)+"','Pa')")



     val sql1:String = "insert into  USA_MODEL values(?,?,?)"
     val preps1 = conn.prepareStatement(sql1)


     for(j in 0..nh-1) {
         val hh: Double = h * j
         for(i in 0..gm.maxgases())
         {

                 preps1.setDouble(1,hh)
                 preps1.setInt(2,i)
                 preps1.setDouble(3,gm.getg(i,0.0,0.0,hh))
                 preps1.addBatch()
         }

         preps1.setDouble(1,hh)
         preps1.setInt(2,-2)
         preps1.setDouble(3,gm.getp(0.0,0.0,hh))
         preps1.addBatch()

         preps1.setDouble(1,hh)
         preps1.setInt(2,-1)
         preps1.setDouble(3,gm.gett(0.0,0.0,hh))
         preps1.addBatch()


     }

     preps1.executeBatch()
     state.close()
     conn.commit()


     conn.close()



 }

}


data class tbasehitran(
    var moln:Int,
    var iso:Int,
    var wnn:Double,
    var Snn:Double,
    var wair:Double,
    var wself:Double,
    var Enn:Double,
    var n:Double,
    var sig:Double,
    var mmol:Double
    )

class funct(cx:Double,cy:Double): function()
{
    var x = 0.0
    var y = 0.0
    init {x = cx
        y = cy
    }
    override  fun func(t:Double):Double
    {
        return exp(-t*t)/(y*y+sqr(x-t))
    }

}

class funct1(cx:Double,cy:Double): function()
{
    var x = 0.0
    var y = 0.0
    init {x = cx
        y = cy
    }
    override  fun func(t:Double):Double
    {
        val l = log(t)
        val ex = exp(-sqr(l))/t
        val sy = sqr(y)

        return ex/(sy+sqr(x+l))

    }

}


class funct2(cx:Double,cy:Double): function()
{
    var x = 0.0
    var y = 0.0
    init {x = cx
        y = cy
    }
    override  fun func(t:Double):Double
    {

        val l = log(t)
        val ex = exp(-sqr(l))/t
        val sy = sqr(y)

        return ex/(sy+sqr(x-l))
    }

}


object Hitran {
    val hitgase: Array<String> = arrayOf("H2O", "CO2", "O3", "N2O", "CO", "CH4", "O2", "NO", "SO2", "NO2", "NH3", "HNO3", "OH", "HF", "HCl", "HBr", "HI", "ClO", "OCS", "H2CO", "HOCl", "N2", "HCN", "CH3Cl", "H2O2", "C2H2", "C2H6", "PH3", "COF2", "SF6", "H2S", "HCOOH", "HO2", "O", "ClONO2", "NO+", "HOBr", "C2H4")
    val coefbolc = 1.380662e-23
    val p1atm = 98066.5
    val aem = 1.6605655e-27
    val Tref = 296.0
    val c2 = 1.4388 //#[cm*K].
    val ermitx: Array<Double> = Array<Double>(10, { 0.0 })
    val ermita: Array<Double> = Array<Double>(10, { 0.0 })
    lateinit var statsum: Array<Array<Double>>
    lateinit var maxt: Array<Double>
    lateinit var mint: Array<Double>
    lateinit var ht: Array<Double>
    lateinit var nt: Array<Int>


    val fxx1 = 1.6
    val fxx2 = 3.55
    val fyy1 = 0.01
    val fyy2 = 0.7
    val fnnx = 100
    val fnny = 2000
    val minfy=0.00000000001
    val fhx1 = (1.6001-0.0)/fnnx
    val fhy1 = (0.7001-0.0199)/fnny

    val fhx2 = (3.5501-1.5999)/fnnx
    val fhy2 = (0.7001-0.0999)/fnny
    val fhx3 = (3.5501-1.5999)/fnnx
    val fhy3 = (0.1-1e-15)/fnny

    val funfxy1: Array<Array<Double>> =  Array(fnny, { Array<Double>(fnnx, {0.0}) })
    val funfxy2: Array<Array<Double>> =  Array(fnny, { Array<Double>(fnnx, {0.0}) })
    val funfxy3: Array<Array<Double>> =  Array(fnny, { Array<Double>(fnnx, {0.0}) })

    val nHit = 38
    val gasesHitM: Array<Double> = Array<Double>(nHit, { 0.0 })

    val Giso: Array<Array<Double>> = arrayOf(
            arrayOf(0.997317, 0.00199983, 0.000372, 0.00031069, 0.000000623, 0.000000116, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.98420, 0.01106, 0.0039471, 0.000734, 0.00004434, 0.00000825, 0.0000039573, 0.00000147, 0.0, 0.0),
            arrayOf(0.992901, 0.00398194, 0.00199097, 0.000740, 0.000370, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.990333, 0.0036409, 0.0036409, 0.00198582, 0.000369, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.98654, 0.01108, 0.0019782, 0.000368, 0.00002222, 0.00000413, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.98827, 0.01110, 0.00061575, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.995262, 0.00399141, 0.000742, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.993974, 0.0036543, 0.00199312, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.94568, 0.04195, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.991616, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.9958715, 0.0036613, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.989110, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.997473, 0.00200014, 0.00015537, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.99984425, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.757587, 0.242257, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.50678, 0.49306, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.99984425, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.75591, 0.24172, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.93739, 0.04158, 0.01053, 0.007399, 0.001880, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.98624, 0.01108, 0.0019776, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.75579, 0.24168, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.9926874, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.98511, 0.01107, 0.0036217, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.74894, 0.23949, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.994952, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.97760, 0.02197, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.97699, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.99953283, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.98654, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.95018, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.94988, 0.04214, 0.007498, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.983898, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.995107, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.997628, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.74957, 0.23970, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.993974, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.5056, 0.4919, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            arrayOf(0.9773, 0.02196, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))


    init {

        val mH = 1.00794
        val mO = 15.9994
        val mC = 12.011
        val mN = 14.00674
        val mF = 18.9984032
        val mP = 30.973762
        val mS = 32.066
        val mCl = 35.4527
        val mBr = 79.904
        val mI = 126.9044
        gasesHitM[0] = aem * mH * 2 + mO * aem
        gasesHitM[1] = aem * mC + aem * mO * 2
        gasesHitM[2] = aem * mO * 3
        gasesHitM[3] = aem * mN * 2 + aem * mO
        gasesHitM[4] = aem * mC + aem * mO
        gasesHitM[5] = aem * mC + aem * mH * 4
        gasesHitM[6] = aem * mO * 2
        gasesHitM[7] = aem * mN + aem * mO
        gasesHitM[8] = aem * mS + aem * mO * 2
        gasesHitM[9] = aem * mN + aem * mO * 2
        gasesHitM[10] = aem * mN + aem * mH * 3
        gasesHitM[11] = aem * mH + aem * mN + aem * mO * 3
        gasesHitM[12] = aem * mO + aem * mH
        gasesHitM[13] = aem * mH + aem * mF
        gasesHitM[14] = aem * mH + aem * mCl
        gasesHitM[15] = aem * mH + aem * mBr
        gasesHitM[16] = aem * mH + aem * mI
        gasesHitM[17] = aem * mCl + aem * mO
        gasesHitM[18] = aem * mO + aem * mC + aem * mS
        gasesHitM[19] = aem * mH * 2 + aem * mC + aem * mO
        gasesHitM[20] = aem * mH + aem * mO + aem * mCl
        gasesHitM[21] = aem * mN * 2
        gasesHitM[22] = aem * mH + aem * mC + aem * mN
        gasesHitM[23] = aem * mC + aem * mH * 3 + aem * mCl
        gasesHitM[24] = aem * mH * 2 + aem * mO * 2
        gasesHitM[25] = aem * mC * 2 + aem * mH * 2
        gasesHitM[26] = aem * mC * 2 + aem * mH * 6
        gasesHitM[27] = aem * mP + aem * mH * 3
        gasesHitM[28] = aem * mC + aem * mO + aem * mF * 2
        gasesHitM[29] = aem * mS + aem * mF * 6
        gasesHitM[30] = aem * mH * 2 + aem * mS
        gasesHitM[31] = aem * mH * 2 + aem * mC + aem * mO * 2
        gasesHitM[32] = aem * mH + aem * mO * 2
        gasesHitM[33] = aem * mO
        gasesHitM[34] = aem * mCl + aem * mO * 3 + aem * mN
        gasesHitM[35] = aem * mN + aem * mO
        gasesHitM[36] = aem * mH + aem * mO + aem * mBr
        gasesHitM[37] = aem * mH * 4 + aem * mC * 2



        val permitx: Array<Double> = arrayOf(3.4361583, 2.5327307, 1.7566830, 1.0366104, 0.3429013, -0.3429013, -1.0366104, -1.7566830, -2.5327307, -3.4361583)
        println("ermita")
        for (i in 0..9)
            ermitx[i] = permitx[i];
        var h0 = 1.0
        for (i in 0..9) {
            h0 = 1.0
            var h1 = 2 * ermitx[i]
            for (j in 2..10) {
                var h2 = 2 * ermitx[i] * h1 - 2 * (j - 1) * h0
                h0 = h1
                h1 = h2
            }
            ermita[i] = (pow(2.0, 9.0)) * fact(9) * sqrt(Math.PI) / (10 * h0 * h0)
            print("  "+ermitx[i]+ " ")
        }


        Class.forName("org.sqlite.JDBC")
        val conn: Connection = DriverManager.getConnection("jdbc:sqlite:common.db")
        val state: Statement = conn.createStatement()
        conn.setAutoCommit(false)





        statsum = Array<Array<Double>>(28, { arrayOf()} )
        nt = Array<Int>(28,{0})
        maxt = Array<Double>(28,{0.0})
        mint = Array<Double>(28,{0.0})
        ht = Array<Double>(28,{0.0})
        for (i in 0..27) {
            var res = state.executeQuery("select count(*) as cnt from STATSUM where num="+i)
            res.next()
            nt[i] = res.getInt("cnt")

            res = state.executeQuery("select max(temp) as maxt from STATSUM where num="+i)
            res.next()
            maxt[i] = res.getDouble("maxt")

            res = state.executeQuery("select min(temp) as mint from STATSUM where num="+i)
            res.next()
            mint[i] = res.getDouble("mint")

            ht[i] = (maxt[i] - mint[i]) / nt[i]
            println(ht[i])
            println(mint[i])
            println(maxt[i])
            println(nt[i])

            statsum[i] = Array<Double>(nt[i],{0.0})

            res = state.executeQuery("select temp, value from statsum where num=" + i + " order by temp")
            var j = 0
            while (res.next() || j < nt[i]) {
                statsum[i][j] = res.getDouble("value")
                j++


            }
        }

    /*for(i in 0..fnny-1) {
        val y = i*fhy1+0.0199
        for (j in 0..fnnx - 1) {
            val x = j*fhx1+0.0
            funfxy1[i][j] = funf1(x, y)

        }

    }

        for(i in 0..fnny-1) {
            val y = i*fhy2+0.0999
            for (j in 0..fnnx - 1) {
                val x = j*fhx2+1.5999
                funfxy2[i][j] = funf1(x, y)

            }

        }*/



    }







    fun readbhit(ar:String, rh:tbasehitran, lam : Array<Double>):Boolean
    {
        val sc = Scanner(ar)
        if(ar.length<160) return false

        rh.moln=sc.findWithinHorizon(".{2}",2).trim().toInt()
        if(rh.moln>26) return false;
        rh.iso=sc.findWithinHorizon(".{1}",1).trim().toInt()
        rh.wnn=sc.findWithinHorizon(".{12}",12).trim().toDouble()
        if(rh.wnn<(lam[0]-100) || rh.wnn>(lam[lam.size-1]+100)) return false


        /*if(ind>=nlam-1) ind = nlam-2;
        if(std::abs(lam[ind]-rh.wnn)<std::abs(lam[ind+1]-rh.wnn))
            if(rh.wnn<(lam[ind]-10) or rh.wnn>(lam[ind]+10)) return false;
        else
            if(rh.wnn<(lam[ind+1]-10) or rh.wnn>(lam[ind+1]+10)) return false;
            */
        rh.Snn=sc.findWithinHorizon(".{10}",10).trim().toDouble()
        sc.findWithinHorizon(".{10}",10)
        rh.wair=sc.findWithinHorizon(".{5}",5).trim().toDouble()
        rh.wself=sc.findWithinHorizon(".{5}",5).trim().toDouble()
        rh.Enn=sc.findWithinHorizon(".{10}",10).trim().toDouble()
        rh.n=sc.findWithinHorizon(".{4}",4).trim().toDouble()
        rh.sig=sc.findWithinHorizon(".{8}",8).trim().toDouble()

        return true;
    }




    fun dopler(bh: tbasehitran, w: Double, pp: Double, T: Double, m: Double): Double {
            val ln2 = sqrt(log(2.0))
            val ldop = bh.wnn * ln2 * sqrt(2 * coefbolc * T / m) / 3e8
            val ss = bh.Snn * exp(-c2 * bh.Enn / T) * (1 - exp(-c2 * bh.wnn / T)) / (exp(-c2 * bh.Enn / Tref) * (1 - exp(-c2 * bh.wnn / Tref)))
            return (p1atm / (T * coefbolc * 1e+6)) * ss * ln2 * exp(-sqr(ln2 * (w - bh.wnn) / ldop)) / (sqrt(Math.PI) * ldop)
        }


    /*val fxx1 = 1.6
    val fxx2 = 3.55
    val fyy1 = 0.01
    val fyy2 = 0.7*/
        fun funf(x: Double, y: Double): Double {
            var val1 = 0.0
            if ((abs(y) >= 0.7) or ((abs(x) >= 3.55) and (abs(y)>=0.01)) or (abs(x)>5.000001) )
                for (i in 0..9)
                    val1 = val1 + ermita[i] / (y * y + sqr(x - ermitx[i]))

            else
                if((abs(y)<=0.02) and  (abs(x)<=1.6)) {
                    val xx = x*x
                    val x4 = xx*xx
                    val yy = y*y
                    val1 = exp(-xx/(y+1))/(yy*0.351+y*0.3183-(x4*1.4-xx*0.38)*yy*0.093)
                   // println(" === "+val1)
                }
                else
                    val1 = funf1(x,y)
                /*if(abs(y)<fyy1)
                    val1 = funfxy1.getxy(abs(x),abs(y),fxx1,fhx1,minfy,fhy1)
                else
                    val1 = funfxy2.getxy(abs(x),abs(y),0.0,fhx2,fyy1,fhy2)*/
                    /*if((abs(y)<0.01) and (abs(x)>=1.6) and (abs(x)<=5.000001))
                    val1 = funf1(x,y)
                    else
                        if((abs(x)<=1.6) and (abs(y)>0.02))
                            val1 = funfxy1.getxy(abs(x),abs(y),0.0,fhx1,0.0199,fhy1)
                        else
                            if((abs(y)>0.1)and(abs(x)>1.6))
                                val1 = funfxy2.getxy(abs(x),abs(y),1.5999,fhx2,0.0999,fhy2)*/


            return val1
        }

        fun funf1(x: Double, y: Double): Double {

            /*val functv1 = funct1(x,y)
            val functv2 = funct2(x,y)

            val a =0.0001
            val b= 0.9999
            return integrate(functv1,a,b,0.01)+integrate(functv2,a,b,0.01)*/
            val funct1 = funct(x,y)
            return integrate(funct1,-5.0,5.0,0.01)

        }

        open fun qt(moln: Int, T: Double): Double {
            if (moln < 28 && moln > 0)
                return statsum[moln - 1].getx(T, mint[moln-1], ht[moln-1])
            else return 1.0
        }

        fun foygt(bh: tbasehitran, w: Double, p: Double, T: Double, ps: Double, m: Double): Double {
            if ((p < 0.00000001)or(T < 1)) return 0.0
            val ln2 = sqrt(log(2.0))
            val ln2pi = sqrt(log(2.0) / Math.PI) / Math.PI
            val ldop = bh.wnn * ln2 * sqrt(2 * coefbolc * T / m) / 3e+8
            val nupt = (pow(Tref / T, bh.n)) * (bh.wair * (p - ps) + bh.wself * ps)
            val ss = (qt(bh.moln, Tref) / qt(bh.moln, T)) * bh.Snn * exp(-c2 * bh.Enn / T) * (1 - exp(-c2 * bh.wnn / T)) / (exp(-c2 * bh.Enn / Tref) * (1 - exp(-c2 * bh.wnn / Tref)))
            val y = nupt * ln2 / ldop

            val foyc = funf((w - bh.wnn) * ln2 / ldop, y) * y * ln2pi / ldop
            return foyc * ss * p1atm / (T * coefbolc * 1.0e+6)
        }

        fun lourence(bh: tbasehitran, www: Double, p: Double, T: Double, ps: Double): Double {
            if ((p < 0.00000001)or(T < 1)) return 0.0

            val nupt = pow(Tref / T, bh.n) * (bh.wair * (p - ps) + bh.wself * ps);
            val ss = (qt(bh.moln, Tref) / qt(bh.moln, T)) * bh.Snn * exp(-c2 * bh.Enn / T) * (1 - exp(-c2 * bh.wnn / T)) / (exp(-c2 * bh.Enn / Tref) * (1 - exp(-c2 * bh.wnn / Tref)))
            val ss1 = ss * nupt * p1atm / (T * coefbolc * 1.0e+6)
            return ss1 / (Math.PI * (sqr(nupt) + sqr(www - (bh.wnn + bh.sig * p))));
        }

        fun contour(bh: tbasehitran, w: Double, pp: Double, T: Double, n: Double, m: Double): Double {
            //double val = foygt(bh,w,pp,T,n,m)/*+foygt(bh,w-0.015,pp,T,n,m)+foygt(bh,w+0.015,pp,T,n,m)*/;
            //return val;
            return foygt(bh, w, pp, T, n, m);
            //return (lourence(bh,w,pp,T,n)+lourence(bh,w-0.015,pp,T,n)+lourence(bh,w+0.015,pp,T,n))/3;
            //return lourence(bh,w,pp,T,n);
        }




    data class tpTc(var p:Double,var T:Double,var c:Array<Double>)
    data class absorpt (var mol:Array<Array<Double>>,var all:Array<Double>)
    //height in km
    fun makeatm(h:Array<Double>):Array<tpTc>
    {
        val pTc:Array<tpTc> = Array(h.size,{tpTc(0.0,0.0,Array(GaseModelDbUsa.maxgases(),{0.0}))})


        for(i in 0..h.size-1) {
            for (j in 0..GaseModelDbUsa.maxgases() - 1)
                pTc[i].c[j] = Ph.Pa_to_atm(GaseModelDbUsa.getg(j, 0.0, 0.0, h[i]))
            pTc[i].p =Ph.Pa_to_atm(GaseModelDbUsa.getp(0.0, 0.0, h[i]))
            pTc[i].T = GaseModelDbUsa.gett(0.0, 0.0, h[i])
        }

        return pTc
    }
    //lam in [cm-1]
    //p in atm
    //T in K
    //c in atm
    fun getallabsorptionpt(lam:Array<Double>,dlam:Double,pTc:Array<tpTc>):Array<absorpt>
    {
        val f:BufferedReader = BufferedReader(FileReader("dat/HITRAN12.par"))
        var s:String = ""
        val rh:tbasehitran= tbasehitran(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
        val a:Array<absorpt> = Array(pTc.size,{absorpt(Array(40,{Array(lam.size,{0.0})}),Array(lam.size,{0.0}))})

        do {
            s=f.readLine()



            if(readbhit(s,rh,lam) )
            {
                //println(""+pTc[0].p+" "+pTc[0].T)

                for(j in 0..pTc.size-1) {



                    for (i in 0..lam.size - 1) {
                        val c = contour(rh, lam[i], pTc[j].p, pTc[j].T, pTc[j].c[rh.moln - 1], gasesHitM[rh.moln - 1])
                        a[j].mol[rh.moln-1][i] +=  c*Giso[rh.moln-1][rh.iso-1]
                        a[j].all[i] +=  c * Giso[rh.moln - 1][rh.iso - 1] * pTc[j].c[rh.moln - 1]

                    }
                }
            }

        } while(s.length>0)

        return a
    }

}



object OpticalModel
{
//p - pressure [atm]; t - temperature [K]; l - wavelength [nm]; result scattering in [km-1]
    fun msct(p:Double,t:Double,l:Double):Double
    {
        val a4 = 1.38e-08
        val t0 = 288.0
        val p0 = 1.013
        if((p<=0.00000000001)or(t<=1)) return 0.0
        val fr1 = l/1000.0
        val pp = p*1013.0
        val cst = (p*t0)/(p0*t)
        val a1=a4*pow((77.6*pp/t)+(0.584*pp/(t*pow(fr1,2))),2)/pow(fr1,4);
        return a1/cst
    }
    //N - concentration [cm-3]; t - temperature [K]; l - wavelength [mkm]; result scattering in [km-1]
    fun msctn(N:Double,t:Double,l:Double):Double = 8.487e+5*pow(Math.PI,3)*pow(pow(Ph.nair(Ph.NinP_air(N,t),t,l),2)-1,2)/(3*N*pow(l*1e-4,4))
    fun msctpin(N:Double,t:Double,l:Double):Double = msctn(N,t,l)*3/(Math.PI*8)

    //#Molecular Backward scatterring [km-1]
    fun msctpi(p:Double,t:Double,l:Double):Double
    {
        val a4 = 1.38e-08
        val t0 = 288.0
        val p0 = 1.013
        val fr1= l/1000.0
        val a3 = 3.0 / (8*Math.PI)
        val pp = p*1013.0
        val cst = (p*t0)/(p0*t)
        val a1=a4*pow((77.6*pp/t)+(0.584*pp/(t*pow(fr1,2))),2)/pow(fr1,4)
        return (a3*a1)/cst
    }



}

class DIAL_IPDA()
{
    var grad:Double = Math.PI
    lateinit var curmolabs:Array<Hitran.absorpt>
    fun makeh(h:Double,r:Double,al:Double,n:Int):Array<Double> = Array(n,{i->h+(r/n)*i*sin(al)})


    fun calc_molabs(lam:Array<Double>,pTc:Array<Hitran.tpTc>):Array<Hitran.absorpt> = Hitran.getallabsorptionpt(lam,0.001,pTc)

    fun calc_aerext(lam:Array<Double>,h:Array<Double>):Array<Array<Double>>
          = Array(lam.size,{i->Array(h.size,{j->AerosolModelIAO.getext(10000.0/lam[i],0.0,0.0,h[j])})})
    fun calc_aersctpi(lam:Array<Double>,h:Array<Double>):Array<Array<Double>>
            = Array(lam.size,{i->Array(h.size,{j->AerosolModelIAO.getsctang(10000.0/lam[i],0.0,0.0,h[j],Math.PI)})})

    fun calc_molsct(lam:Array<Double>,pTc:Array<Hitran.tpTc>):Array<Array<Double>> =
            Array(lam.size,{i->Array(pTc.size,{j->OpticalModel.msct(pTc[j].p,pTc[j].T,1.0e+7/lam[i])})})
    fun calc_molsctpi(lam:Array<Double>,pTc:Array<Hitran.tpTc>):Array<Array<Double>> =
            Array(lam.size,{i->Array(pTc.size,{j->OpticalModel.msctpi(pTc[j].p,pTc[j].T,1.0e+7/lam[i])})})
    data class signal(var sig:Double,var tau:Double,var molabs:Double,var bpi:Double)

    fun norm_noise():Double
    {
        var r:Double = 0.0
        for(i in 0..5)
            r = r+random()
        r = (r - 3)/3
        return r

    }

    fun noise_mean(p:Double,amp:Double):Double
    {
        val v = norm_noise()
        var n:Double
        if((v<0) && (p<amp))  n = p+v*p
        else n = p+v*amp
        return n
    }





    fun calc_snr()
    {
        val E0=100.0e-3
        val r0 = 0.5e-3
        val A = Math.PI*r0*r0
        val nu = 0.8
        val stdist = 450.0
        val M = 20.0
        val Q = 0.8
        val B = 20e+6
        val Rr = 1.014
        val F = 4.3
        val id = 160e-15
        val i0 = 4e-15
        val u0 = 3e-9
        val T = 290
        val Rf = 1e+6
        val  Cd = 4e-12
        //val Pb = 1.2e-9
        val e = 1.602176565e-19
        val kb = 1.3806488e-23
        val hn = 6.62606E-34
        val nu1 = 3e+8/1.5e-6
        val Nimp = 10000
        val c=3.0e+8/1.0e+3
        val dh = 0.1
        val sigon = 1



        val H=23*1000
        val ro = 0.31
        val E = 0.25 //W/m2*nm
        val dl = 1

        val ug = 10.0*1e-3
        val rp = H*sin(ug)/cos(ug)
        val Sp = sqr(rp)*Math.PI
        var Pb = E*ro*dl*A*Sp/sqr(H)

        println("Pb= "+ Pb)



        val n = (H*1e-3/dh).toInt()
        val lam = arrayOf(1e+7/1572.018,1e+7/1572.19)
        val h = makeh(H*1e-3,H*1e-3,-Math.PI/2,n)
        val sig =   calc_signal(1,lam,h,dh)

        var Pb1 = 0.0
        var Pb2 = 0.0
        for(i in 2..n-1) {
            val rp = h[i]*tan(ug)
            val Sp = sqr(rp)*Math.PI
            val Pb11 = E*sig[0][i-1].bpi*dh*dl*Sp*A/sqr(h[i-1])
            val Pb12 = E*sig[0][i].bpi*dh*dl*Sp*A/sqr(h[i])
            Pb1 = Pb1+(Pb11+Pb12)*0.5
            val Pb21 = E*sig[1][i-1].bpi*dh*dl*Sp*A/sqr(h[i-1])
            val Pb22 = E*sig[1][i].bpi*dh*dl*Sp*A/sqr(h[i])
            Pb2 = Pb2+(Pb21+Pb22)*0.5


        }

        val Pba = (Pb1+Pb2)/2

        println("Pba = "+ Pba)

        Pb = Pb + Pba

        val Pon = Array<Double>(n,{0.0})
        val Poff = Array<Double>(n,{0.0})
        val cnron = Array<Double>(n,{0.0})
        val cnroff = Array<Double>(n,{0.0})
        val noise = Array<Double>(n,{0.0})


        for(i in 0..n-1) {
            Pon[i] = sig[0][i].sig * dh * A * Q * nu * E0/(dh/c)
            Poff[i] = sig[1][i].sig * dh * A * Q * nu * E0/(dh/c)

            val P1 = Pon[i]
            val P2 = Poff[i]
            val i2non = B*(2*e*M*M*F*Rr*(P1+Pb)+id*id+i0*i0+4*kb*T/Rf+pow(u0/Rf,2))+pow(B,3)/3*pow(u0*2*Math.PI*Cd,2)
            val i2noff = B*(2*e*M*M*F*Rr*(P2+Pb)+id*id+i0*i0+4*kb*T/Rf+pow(u0/Rf,2))+pow(B,3)/3*pow(u0*2*Math.PI*Cd,2)
            cnron[i] = (P1*M*Rr/sqrt(i2non))*sqrt(Nimp.toDouble())
            cnroff[i] = (P2*M*Rr/sqrt(i2noff))*sqrt(Nimp.toDouble())
            noise[i] = (Pon[i]/cnron[i]+Poff[i]/cnroff[i])/2
            println(" " +Pon[i]+" "+noise[i])
        }

        val conc = Array<Double>(n,{0.0})
        val conc1 = Array<Double>(n,{0.0})
        val difc = Array<Double>(n,{0.0})
        val err = Array<Double>(n,{0.0})
        val p = Hitran.makeatm(h)

        val Nerr = 10000
        for(k in 0..Nerr-1) {
            if(k % 1000 == 0) println(k)
            val Pone = Array<Double>(n,{0.0})
            val Poffe = Array<Double>(n,{0.0})
            val Ni = 1

            for (i in 0..n - 1) {
                for (j in 0..Ni - 1) {
                    Pone[i] += noise_mean(Pon[i], noise[i])
                    Poffe[i] += noise_mean(Poff[i], noise[i])

                }
                Pone[i] = Pone[i] / Ni
                Poffe[i] = Poffe[i] / Ni
            }

            for (i in 0..n - 1) {
                difc[i] = -log(Pone[i] / Poffe[i])
            }

            val sdif1 = difc.smootmean(8)
            val sdif = sdif1.smootmean(8)


            for (i in 0..n - 2) {
                conc[i] = (sdif[i + 1] - sdif[i]) / ((sig[0][i].molabs - sig[1][i].molabs) * dh * 2 * 1e+5)
                conc1[i] = p[i].c[1]
                err[i] += abs(conc[i] - conc1[i]) / conc1[i]

            }

        }

        for (i in 0..n - 2)
            err[i] = err[i]/Nerr


        val f = FileWriter("dat1/f.dat")

        val err1 = err.smootmean(1)
        val err2 = err1.smootmean(1)
        for(i in 0..n-1)
            f.write(""+h[i]+"  "+err2[i]*100.0+"\r\n")

        f.close()





    }

    fun calc_signal(nmol:Int,lam:Array<Double>,h:Array<Double>,dr:Double):Array<Array<signal>>
    {
        val sig:Array<Array<signal>> = Array(lam.size,{Array(h.size,{signal(0.0,0.0,0.0,0.0)})})

        val pTc = Hitran.makeatm(h)
        val absr = calc_molabs(lam,pTc)
        val aerext = calc_aerext(lam,h)
        val aerpi = calc_aersctpi(lam,h)
        val msct = calc_molsct(lam,pTc)
        val msctpi = calc_molsctpi(lam,pTc)

        for(i in 0..lam.size-1)
        {
            for(j in 1..h.size-1) {
                sig[i][j].tau = (absr[j-1].all[i] + absr[j].all[i])*1.0e+5*0.5*dr
                sig[i][j].tau += (aerext[i][j-1]+aerext[i][j])*0.5*dr
                sig[i][j].tau += (msct[i][j-1]+msct[i][j])*0.5*dr
                sig[i][j].tau =sig[i][j-1].tau+ sig[i][j].tau
                sig[i][j].bpi = (aerpi[i][j]+msctpi[i][j])
                sig[i][j].sig = (sig[i][j].bpi) * exp(-2.0*sig[i][j].tau)/sqr(dr*j)
                sig[i][j].molabs = absr[j].mol[nmol][i]
            }
            sig[i][0].sig = sig[i][1].sig*sqr(dr)/sqr(dr/3)
            sig[i][0].molabs = absr[0].mol[nmol][i]
            sig[i][0].bpi = (aerpi[i][0]+msctpi[i][0])
        }


        return sig
    }

    fun calcP_IPDA(lam:Double):Double
    {
        var P:Double = 0.0

        return P
    }

    fun calcP_DIAL():Array<Double>
    {

        return arrayOf(0.0)
    }

}