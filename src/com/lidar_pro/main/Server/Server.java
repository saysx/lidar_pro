package com.lidar_pro.main.Server;

import com.lidar_pro.main.calc.LidarKt;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;
import java.util.ArrayList;

import static com.lidar_pro.main.Helpers.GlobalVars.*;

/**
 * Created by Iwan on 27.07.2016.
 */
public class Server implements Runnable {

    Socket clientSocket;

    BufferedReader in = null;
    PrintWriter out = null;

    String input = null;
    String output = null;

    public Server(Socket clientSocket) {
        this.clientSocket = clientSocket;
    }

    @Override
    public void run() {
        try {
            in = new BufferedReader(new InputStreamReader(clientSocket.getInputStream())); // read from client
            out = new PrintWriter(clientSocket.getOutputStream(), true);// send response to client

            while((input = in.readLine()) != null) {
                //out.println("S ::: " + input);
                //System.out.println(input);
                String s[] = input.split(","); // потому что в аррейлист нельзя разбить стринг по делиметрам
                ArrayList<String> al = new ArrayList<String>(); // потому что можно проверить, существует ли тот или иной элемент в массиве

                for(int i = 0; i < s.length; i++) {
                    al.add(s[i]);
                }
                for (String ss : al) {
                    System.out.println("Incoming request: " + ss);
                }
                if (al.get(0).equals("abs")) {
                    LidarKt.calc_absorption_spectre(Double.parseDouble(al.get(1)), Integer.parseInt(al.get(2)));
                    out.println(abs_fileName + "," + abs_ms_fileName);
                }

                if (al.get(0).equals("snr")) {
                    //if (al.size() > 3) { // smootmean
                    System.out.println(al.size());
                    for (int i = 3; i < al.size(); i++ ) {
                        if (al.get(i).length() > 4) {
                            if (al.get(i).substring(0, 4).equals(smoothing_prefix)) {
                                smootMean_coef = Integer.parseInt(al.get(i).substring(4, 5));
                            }
                            if (al.get(i).substring(0, 4).equals(calibration_coef_prefix)) {
                                calibration_coef = (al.get(i).substring(4, 5));
                            }
                        }
                    }
                    LidarKt.calc_dial_snr(Double.parseDouble(al.get(1)), Double.parseDouble(al.get(2)), Integer.parseInt(al.get(3)), smootMean_coef);
                    out.println(snr_fileName);
                }
            }
            /*out.close();
            in.close();
            clientSocket.close();*/
        } catch(IOException e) {
            e.printStackTrace();
        }
    }
}
