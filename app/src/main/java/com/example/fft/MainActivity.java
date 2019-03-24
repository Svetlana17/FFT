package com.example.fft;

import android.annotation.SuppressLint;
import android.content.Context;
import android.content.Intent;
import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;
import android.net.Uri;
import android.os.Build;
import android.support.v7.app.AppCompatActivity;
import android.os.Bundle;
import android.text.Editable;
import android.text.TextWatcher;
import android.util.Log;
import android.view.MotionEvent;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;
import android.widget.Toast;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.FloatBuffer;
import java.sql.Time;
import java.util.LinkedList;
import java.util.Timer;
import java.util.TimerTask;
import java.util.concurrent.TimeUnit;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import static android.util.Half.EPSILON;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

public class MainActivity extends AppCompatActivity  implements SensorEventListener {



        boolean isRunning;
        final String TAG = "SensorLog";
        EditText editTextShag;
        int v;
        private FloatBuffer[] inFloatBuffer;
        private FFT fft;
        private float frequency;
        private SensorManager sensorManager;
        private Sensor accelerometer;
        private float[] acceleration;
        private float magnitude;
        private float tMean;
        private float tStd;
        private float time;
        private float timeLast;
        private float timeClosestOnset;
        private float timeNextOnset;
        private float timeOnsetLast;
        private boolean tracking;
        private boolean onset;
        private boolean onsetFired;
        private static final int NUM_CHANNEL_PLUS_ONE = 4;
        private int fftSize;
        private int fftHopSize;
        private int fftSizePadded;
        private int fftSizePaddedCompact;
        private static final float FFT_TIME = 4f;
        private long timeMs;
        private long timeLastMs;
        private long timeStartMs;
        private float sumN;
        private float sum;
        private float sumOfSquared;
        private float ratio;
        private static final float MILLIS_TO_SECONDS = 1 / (float) TimeUnit.SECONDS.toMillis(1);
        private static final int NUM_CHANNEL = 3;
        private static final int FFT_CUTOFF_BPM = 60;
        private static final float RATIO_CUTOFF = 20;
        @SuppressLint("ClickableViewAccessibility")
        @Override
        protected void onCreate(Bundle savedInstanceState) {
            super.onCreate(savedInstanceState);
            //  setContentView(R.layout.activity_mainp);
            setContentView(R.layout.activity_main);

            editTextShag.addTextChangedListener(new TextWatcher() {
                @Override
                public void beforeTextChanged(CharSequence charSequence, int i, int i1, int i2) {
                }

                @Override
                public void onTextChanged(CharSequence charSequence, int i, int i1, int i2) {
                }

                @Override
                public void afterTextChanged(Editable editable) {
                    v = Integer.parseInt(editable.toString());
                }
            });

            isRunning = false;

            sensorManager = (SensorManager) getSystemService(Context.SENSOR_SERVICE);



///

        }

    public MainActivity(boolean isRunning) {
        this.isRunning = isRunning;
    }

    SensorEventListener listener = new SensorEventListener() {


            @Override
            public void onSensorChanged(SensorEvent hardEvent) {
        if(hardEvent.sensor.getType() != Sensor.TYPE_ACCELEROMETER)
        return;

        timeMs = System.currentTimeMillis();
        acceleration = hardEvent.values;

        processTime();
        processOnset();

        processAcceleration();
            }

        @Override
        public void onAccuracyChanged(Sensor sensor, int accuracy) {

        }

        private void processOnset() {
            onset = !onsetFired && (timeClosestOnset != 0) && (timeLast < timeOnsetLast) && (time > timeClosestOnset);
            onsetFired |= onset;
            timeOnsetLast = timeClosestOnset;

        }
        public float getRatio() {
            return ratio;
        }

        private void processAcceleration() {
            magnitude = 0;
            for (int i = 0; i < NUM_CHANNEL; i++) {
                magnitude += (acceleration[i] * acceleration[i]);
                inFloatBuffer[i].put(acceleration[i]);
            }
            inFloatBuffer[NUM_CHANNEL].put(magnitude);

            if (inFloatBuffer[0].position() >= fftSize) {
                float[] inReal = null;
                float[] inImag = null;
                float[] fftSquaredSum = new float[fftSizePaddedCompact];

                for (int i = 0; i < NUM_CHANNEL_PLUS_ONE; i++) {
                    inReal = getFrame(inFloatBuffer[i]);
                    inImag = new float[fftSizePadded];

                    Util.center(inReal, fftSize);
                    fft.fft(inReal, inImag, true);

                    if (i != NUM_CHANNEL)
                        for (int j = 0; j < fftSizePaddedCompact; j++)
                            fftSquaredSum[j] += (inReal[j] * inReal[j] + inImag[j] * inImag[j]);
                }

                int cutOffIndex = (int) Math.ceil(Util.bpmToFreq(FFT_CUTOFF_BPM) * tMean * fftSizePadded);
                ;
                float[] maxStat = Util.max(fftSquaredSum, cutOffIndex);
                int maxIndex = (int) maxStat[1];

                frequency = maxStat[1] / fftSizePadded / tMean;
                ratio = maxStat[0] / Util.mean(fftSquaredSum);
                tracking = (maxIndex > cutOffIndex) && (ratio > RATIO_CUTOFF);

                if (tracking) {
                    float maxAngularFrequency = (float) (2 * PI * frequency);
                    float maxPhase = (float) Math.atan2(inImag[maxIndex], inReal[maxIndex]);
                    float timeToNextBeatFromShifted = -Util.smallestMagnitudeCoterminalAngle(maxPhase + maxAngularFrequency * FFT_TIME) / maxAngularFrequency;

                    timeClosestOnset = time + timeToNextBeatFromShifted;
                    timeNextOnset = timeClosestOnset + ((timeToNextBeatFromShifted > 0) ? 0 : 1 / frequency);
                    onsetFired = false;
                }





            }
        }


        };

    private float[] getFrame(FloatBuffer buffer) {

            float[] in = new float[fftSizePadded];

            buffer.flip();
            buffer.get(in, 0, fftSize);
            buffer.position(fftHopSize);
            buffer.compact();

            return in;
        }


    private void processTime() {

            time = (timeMs - timeStartMs) * MILLIS_TO_SECONDS;
            timeLast = (timeLastMs - timeStartMs) * MILLIS_TO_SECONDS;

            if(timeLastMs != 0) {
                float diffTime = (timeMs - timeLastMs) * MILLIS_TO_SECONDS;

                sumN += 1;
                sum += diffTime;
                sumOfSquared += diffTime * diffTime;

                tMean = sum / sumN;
                tStd = (float) Math.sqrt(sumOfSquared / sumN - tMean * tMean);
            }
            timeLastMs = timeMs;
        }




    @Override
    public void onSensorChanged(SensorEvent event) {

    }

    @Override
    public void onAccuracyChanged(Sensor sensor, int accuracy) {

    }
}
