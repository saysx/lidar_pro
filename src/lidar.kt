/**
 * Created by say on 06.01.16.
 */
package demo

import mie.DIAL_IPDA

data class vas(val v1:Double,val v2:Double)
{

}


fun ret():Array<vas>
{

    return arrayOf(vas(2.0,2.0))
}

fun main(args : Array<String>) {
    val today = java.lang.System.currentTimeMillis()



    /*val pTc:Array<Hitran.tpTc> = Hitran.makeatm(arrayOf(0.0,1.0,2.0))

    println("Pressure of gases - ")
    pTc[0].c[0] =pTc[0].c[0]*2.0
    for(i in 0..26)
        println(pTc[0].c[i]*100)

    val lam = Array<Double>(1000,{i->10000.0+i.toDouble()/100.0})
    val c = Hitran.getallabsorptionpt(lam,1,0.001,pTc)
    val f = FileWriter("dat1/f.dat")

    println(c[0].all.size)
    for(i in 0..c[0].all.size-1)
      f.write(""+lam[i]+"  "+c[0].all[i]+"\r\n")

    f.close()*/

    val di:DIAL_IPDA = DIAL_IPDA()
    /*val h = di.makeh(0.0,100.0,Math.PI/2,200)
    println("height == == ="+h[0])
    println("height == == ="+h[1])
    println(h[199])
    val lam = Array(200,{i->10000.0*1e+3/(1571.8+i*0.5/200)})
    val sig = di.calc_signal(1,lam,h,0.5)

    val f = FileWriter("dat1/f.dat")

    var srex = 0.0
    for(i in 0..sig.size-1) {

        val ex = exp(-sig[i][sig[i].size - 1].tau)
        srex = srex + ex
        f.write("" + lam[i] + "  " + ex + "\r\n")

    }

    srex = srex / sig.size


    f.close()
    println(srex)*/

    /*val r = 1e-3
    val A = sqr(r)*Math.PI
    val H=23*1000
    val ro = 0.31
    val E = 0.25 //W/m2*nm
    val dl = 1

    val ug = 10*1e-3
    val rp = H*sin(ug)/cos(ug)
    val Sp = sqr(rp)*Math.PI
    val Pb = E*ro*dl*A*Sp/sqr(H)
    println(rp)
    println(Sp)
    println(Pb)*/

    di.calc_snr()


    //Hitran.getabsorption(arrayOf(10000.0,10010.0),1,0.1)
   /* val b: blocks = blocks()
    val wr2: writedata2d = writedata2d()
    val str: stringconst = stringconst("dat1/FileName0.dat")
    val intv: intconst = intconst(200)
    val doublev1: doubleconst = doubleconst(3000.0)
    val doublev2: doubleconst = doubleconst(4000.0)
    val doublev3: doubleconst = doubleconst(0.0)
    val aerextv : aerextvec = aerextvec()
    val makev:makevec = makevec()
    val concv:concatvectors = concatvectors()

    b.addblock(wr2)
    b.addblock(str)
    b.addblock(intv)
    b.addblock(doublev1)
    b.addblock(doublev2)
    b.addblock(doublev3)
    b.addblock(aerextv)
    b.addblock(makev)
    b.addblock(concv)

    b.connect(wr2,1,str,0)
    b.connect(wr2,0,concv,0)
    b.connect(concv,0,makev,0)

    b.connect(concv,1,aerextv,0)
    b.connect(makev,0,doublev1,0)
    b.connect(makev,1,doublev2,0)
    b.connect(makev,2,intv,0)
    b.connect(aerextv,0,doublev3,0)
    b.connect(aerextv,1,makev,0)





    b.listb[0].calculate()

*/
    /*val act:logisticact = logisticact()
    val nt: backppg =  backppg(2,arrayOf(10,10,1),arrayOf(act,act,act))
    val n = 100
    val xx:Array<Array<Double>> = Array(n,{arrayOf(random()-0.5,random()-0.5)})
    val yy:Array<Array<Double>> = Array(n,{i->arrayOf(xx[i][0]+xx[i][1])})
    nt.al = 0.5

    for(j in 0..10000)
    for(i in 0..n-1)
    {
        nt.learn(xx[i],yy[i])

    }

    var y = nt.ask(xx[1])

    println(y[0])
    println(yy[1][0])

    y = nt.ask(xx[3])

    println(y[0])
    println(yy[3][0])*/
    /*val lgst:logisticact = logisticact()
    val ne = 50
    val inv: inversetrain = inversetrain(3,20,ne,lgst)
    for(i in 0..ne-1)
    {
     inv.x[0][i] = random()
     inv.x[1][i] = random()
     inv.x[2][i] =inv.x[1][i]*inv.x[0][i]
    }

   // for(i in 0..60)
       // inv.learnstep()


    inv.gradientbatch(0.5,10000)

    println(inv.dataseterror())

    val xx= inv.askxyx(getvector(inv.x,0))
    for(i in 0..inv.x.size-1)
        println(" " + inv.x[i][0]+ "   "+xx[i] )
        */

  /*  val train: trainset = trainset(100,3,3)
    for(i in 0..train.ntrain-1)
    {
        train.set[i][0][0] = random()
        train.set[i][0][1] = random()
        train.set[i][0][2] = train.set[i][0][0]*train.set[i][0][1]

        train.set[i][1][0] = train.set[i][0][0]
        train.set[i][1][1] = train.set[i][0][1]
        train.set[i][1][2] = train.set[i][0][2]

    }
    val lgst:logisticact = logisticact()
    val pl: prelearngradient = prelearngradient(3, arrayOf(6, 6, 6, 6, 3), arrayOf(lgst, lgst, lgst, lgst, lgst), train)
    //val pl: backppg = backppg(3, arrayOf(15, 15, 15, 5, 3), arrayOf(lgst, lgst, lgst, lgst, lgst))
    pl.al = 0.01
    for(j in 0..10000)
        for(i in 0..train.ntrain-1)
        {
            val j = (random()*train.ntrain).toInt()
            pl.learn(train.set[j][0],train.set[j][1])

        }

    for( i in 0..10) {
        var y = pl.ask(train.set[i][0])
        println("x1 = " + y[0] + "  " + train.set[i][1][0])
        println("x2 = " + y[1] + "  " + train.set[i][1][1])
        println("y = " + y[2] + "  " + train.set[i][1][2])
    }*/


    val tm = System.currentTimeMillis() - today
        println(tm)


    }