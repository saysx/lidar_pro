package calculate

import maths.random
import mie.AerosolModelIAO
import java.io.FileWriter

/**
 * Created by say on 19.03.16.
 */
class calculate
{


}

open  class data(cn:String)
{
    open  var name:String = cn
    open  var s:Int = 0
    open  var type:String = ""
    open  var state:Boolean = false





    fun isconnect(x:data):Boolean = when(x)
    {
        is data1d -> when(this) {is data1d -> true else -> false}
        is data2d -> when(this) {is data2d -> true else -> false}
        is data0d -> when(this) {is data0d -> true else -> false}
        is data0s -> when(this) {is data0s -> true else -> false}
        is data0i -> when(this) {is data0i -> true else -> false}
        else -> false

    }





}


open class datad(cn:String):data(cn)
{
    open  var em:String = ""
    init { type = "Double" }
    fun checkem(x:datad) =
            if(x.em == this.em) true else false

}

open class datai(cn:String):data(cn)
{
    open  var em:String = ""
    init { type = "Integer" }
    fun checkem(x:datai) =
            if(x.em == this.em) true else false

}


open class datas(cn:String):data(cn)
{
    init { type = "String" }
}


class data0d(cn:String):datad(cn)
{
    var d:Double = 0.0
    init {
        s = 0

    }

}

class data1d(cn:String):datad(cn)
{
    var d:Array<Double> = Array<Double>(0,{0.0})
    init {
        s = 1

    }

}

class data2d(cn:String):datad(cn)
{
    var d:Array<Array<Double>> = Array<Array<Double>>(0,{Array<Double>(0,{0.0})})
    init {
        s = 2
    }

}

class data0s(cn:String):datas(cn)
{
    var d:String = ""
    init {
        s = 0

    }

}

class data0i(cn:String):datai(cn)
{
    var d:Int = 0
    init {
        s = 0

    }

}


open class block(cdin:Array<data>,cdout:Array<data>)
{

    var din:Array<data> = cdin
    var pdin:Array<data> = Array<data>(cdin.size,{i->data0s("")})
    var pbl:Array<block> = Array<block>(cdin.size,{i->block(arrayOf(),arrayOf())})
    var dout:Array<data> = cdout
    var nin = cdin.size
    var nout = cdout.size
    var allconnect:Boolean = false
    var connections:Array<Boolean> = Array<Boolean>(nin,{false})
    var active = true
    var index:Int = 0
    var name:String = "block"
    init {if (cdin.size == 0) allconnect = true}
    fun disconnect(cnin:Int):Boolean
    {
        if(cnin<=0) return false
        if(cnin>=nin) return false
        if(connections[cnin]) {connections[cnin] = false
            allconnect = false } else return false

        return true
    }

    fun connect(cnin:Int,cnout:Int,tb:block):Boolean
    {
        if(cnin<0) return false
        if(cnin>=nin) return false
        if(cnout<0) return false
        if(cnout>=tb.nout) return false


        if(din[cnin].isconnect(tb.dout[cnout]))
        {
            pbl[cnin]  = tb
            pdin[cnin] = tb.dout[cnout]
            connections[cnin] = true


        } else return false

        var n = 0
        for(i in 0..nin-1)
         {if(connections[i]) n=n+1}
        if(n>=cnin) allconnect = true else allconnect = false

        return true

    }
    open fun childexecute():Boolean = true

    fun execute():Boolean {
        if(allconnect)
        {
            println(name+" execute")
            for(i in 0..nout-1)
                dout[i].state = true
            return childexecute()

        } else return false

    }

    fun calculate():Boolean
    {

        if(active)
        for(i in 0..pbl.size-1)
            //if(pdin[i].state)
                pbl[i].calculate()


        return execute()
    }

}

class writedata1d: block
{

    constructor() : super(arrayOf(data1d(""),data0s("FileName")),arrayOf())
    {name = "writedata1d"}
    override fun childexecute():Boolean {


        val val1 = true
        val in1:data1d = pdin[0] as data1d
        val in2:data0s = pdin[1] as data0s
        val f = FileWriter(in2.d)




        for(i in 0..in1.d.size-1) {
            f.write("" + in1.d[i] + "\r\n")
        }
        f.close()
        return val1
    }

}





class writedata2d: block
{

    constructor() : super(arrayOf(data2d(""),data0s("FileName")),arrayOf())
    {name = "writedata2d"}
    override fun childexecute():Boolean {


        val val1 = true
        val in1:data2d = pdin[0] as data2d
        val in2:data0s = pdin[1] as data0s
        val f = FileWriter(in2.d)

        var min:Int = in1.d[0].size
        println(min)
        println(in1.d[1].size)
        for(i in 1..in1.d.size-1)
            if(in1.d[i].size<min) min =in1.d[i].size

        println(min)

        for(i in 0..min-1) {
            for(j in 0..in1.d.size-1)
             f.write(""+in1.d[j][i]+"\t")
            f.write("\r\n")
            }
        f.close()
        return val1
    }

}






class randomvecd: block
{
    constructor() : super(arrayOf(data0i("VectorSize")),arrayOf(data1d("RandomVector")))
    {name = "randomvecd"}


    override fun childexecute():Boolean
    {
        val val1 = true
        val in1:data0i  = pdin[0] as data0i
        val o1:data1d  = dout[0] as data1d
        println("randsize =  "+in1.d)
        o1.d = Array<Double>(in1.d,{i->random()})
        return val1

    }

}


class intconst: block
{
    var value:Int = 0
    constructor(cval:Int) : super(arrayOf(),arrayOf(data0i("IntegerConst")))
    {
        value = cval
        name = "intconst"

    }
    override fun childexecute():Boolean {
        println(name)
        val o1:data0i  = dout[0] as data0i
        o1.d = value

        return true
    }

}

class doubleconst: block
{
    var value:Double = 0.0
    constructor(cval:Double) : super(arrayOf(),arrayOf(data0d("DoubleConst")))
    {
        value = cval
        name = "doubleconst"

    }
    override fun childexecute():Boolean {
        println(name)
        val o1:data0d  = dout[0] as data0d
        o1.d = value

        return true
    }

}



class stringconst: block
{
    var value:String = ""
    constructor(cval:String) : super(arrayOf(),arrayOf(data0s("StringConst")))
    {
        value = cval
        name = "stringconst"
    }
    override fun childexecute():Boolean {
        val o1:data0s  = dout[0] as data0s
        o1.d = value
        return true
    }


}


class concatvectors:block
{
    constructor() : super(arrayOf(data1d("vector 1"),
            data1d("vector 2")),
            arrayOf(data2d("DoubleMatrix")))
    {
        name = "concatvectors"

    }
    override fun childexecute():Boolean {
        val o1:data2d  = dout[0] as data2d
        val i1:data1d  = pdin[0] as data1d
        val i2:data1d  = pdin[1] as data1d


        o1.d=arrayOf()
        o1.d = o1.d+i1.d+i2.d
        return true
    }


}

class concatvecmatr:block
{
    constructor() : super(arrayOf(data2d("matrix 1"),
            data1d("vector 1")),
            arrayOf(data2d("DoubleMatrix")))
    {
        name = "concatvecmatr"

    }
    override fun childexecute():Boolean {
        val o1:data2d  = dout[0] as data2d
        val i1:data2d  = pdin[0] as data2d
        val i2:data1d  = pdin[1] as data1d


        o1.d = i1.d+i2.d
        return true
    }


}



class makevec:block
{
    constructor() : super(arrayOf(data0d("a border"),
            data0d("b border"),data0i("Size of vector")),
            arrayOf(data1d("DoubleVector")))
    {
        name = "aerextvec"

    }
    override fun childexecute():Boolean {
        val o1:data1d  = dout[0] as data1d
        val i1:data0d  = pdin[0] as data0d
        val i2:data0d  = pdin[1] as data0d
        val i3:data0i  = pdin[2] as data0i
        val h = (i2.d-i1.d)/i3.d
        o1.d = Array<Double>(i3.d,{i->i1.d+h*i})
        return true
    }


}







class aerextvec:block
{
    constructor() : super(arrayOf(data0d("Height in [km]"),data1d("Wavelength Vector in [mkm]")),
                     arrayOf(data1d("Aerosol Extinction in [km-1]")))
    {
       name = "aerextvec"

    }
    override fun childexecute():Boolean {
        val o1:data1d  = dout[0] as data1d
        val i1:data0d  = pdin[0] as data0d
        val i2:data1d  = pdin[1] as data1d


        println(i2.d.size)
        o1.d = Array<Double>(i2.d.size,{i->AerosolModelIAO.getext(i2.d[i],0.0,0.0,i1.d)})

        return true
    }


}


class blocks()
{
    var listb:Array<block> = arrayOf()
    fun addblock(bl:block)
    {
        listb = listb+bl
        listb[listb.size-1].index = listb.size-1
    }
    fun connect(i1:Int,cin:Int,i2:Int,cout:Int):Boolean
    {
        if((i1<0) or (i2<0) or (i1>=listb.size) or (i2 > listb.size)) return false
        return listb[i1].connect(cin,cout,listb[i2])

    }

    fun connect(blin:block,cin:Int,blout:block,cout:Int):Boolean
    {
        return blin.connect(cin,cout,blout)

    }


    fun calculate():Boolean
    {
        return true
    }



}