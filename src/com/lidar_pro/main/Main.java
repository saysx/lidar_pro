package com.lidar_pro.main;

import com.lidar_pro.main.Server.Server;

import java.net.ServerSocket;
import java.net.Socket;

/**
 * Created by Iwan on 27.07.2016.
 */
public class Main {

    /** LWCH Server work algorithm:
     *
     * client's request for generate abs/snr file -> check if file exists by file name in archive:
     * if exists -> send to client file name
     * else -> create file & generate data & suspend work with this file -> send to client file name
     *
     */


    public Main() {
        // TODO: Administrator's UI for monitoring:
        // TODO: User activity log system
        // ?TODO: Archive of generated files

        // Multithreaded server:
        try {
            ServerSocket serverSocket = new ServerSocket(9096);
            while (true) {
                Socket sock = serverSocket.accept();
                System.out.println("Client connected from: " + sock.getInetAddress() + ":" + sock.getPort());
                new Thread(new Server(sock)).start();
            }
        } catch(Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String [  ] laser) {
        new Main();
    }

}




