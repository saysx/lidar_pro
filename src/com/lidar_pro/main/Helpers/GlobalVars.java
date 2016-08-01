package com.lidar_pro.main.Helpers;

/**
 * Created by Iwan on 28.07.2016.
 */
public class GlobalVars {

    public static String abs_fileName = null;
    public static String snr_fileName = null;
    public static String abs_ms_fileName = null;

    public static String spectreMols_text = null;


    public static String fileExtension_dat = ".dat";

    // plugin variables:
    public static String abs_plugin_suffix = "";
    public static String snr_plugin_suffix = "";

    public static int smootMean_coef = 0;

    // plugin variables:

    public static String smoothing_coef = "";
    public static String calibration_coef = "";

    // plugin prefix codes:
// TODO: перенести префиксы в загружаемый файл для парса сервером и клиентом, чтобы не синхронизировать постоянно в коде
    public static String smoothing_prefix = "smth";
    public static String calibration_coef_prefix = "cali";
}
