package com.example.fft;

//package com.example.fft;
//
public class FFT
{
//	/* FIELDS ******************************************************************/
//    private int n;
//    private int m;
//    private int k;
//
//	// The results of the FFT. // Результаты БПФ.
//	private	double[]	real_output;
//	private	double[]	imaginary_output; //мнимая часть
//
//	// The phase angles Фазовые углы
//	private double[]	output_angle;
//
//	// Magnitude and power spectra // Амплитуда и спектры мощности
//	private double[]	output_magnitude;
//	private double[]	output_power;
//
//
//	/* CONSTRUCTOR *************************************************************/
//
//
//	/**
//	 * Performs the Fourier transform and stores the real and imaginary results.
//	 * Input signals are zero-padded if they do not have a length equal to a
//	 * power of 2.
//	 *
//	 * @param	real_input			The real part of the signal to be transformed.
//	 * @param	imaginary_input		The imaginary part of the signal to be.
//	 *								transformed. This may be null if the signal
//	 *								is entirely real.
//	 * @param	inverse_transform	A value of false implies that a forward
//	 *								transform is to be applied, and a value of
//	 *								true means that an inverse transform is tob
//	 *								be applied.
//	 * @param	use_hanning_window	A value of true means that a Hanning window
//	 *								will be applied to the real_input. A value
//	 *								of valse will result in the application of
//	 *								a Hanning window.
//	 * @throws	Exception			Throws an exception if the real and imaginary
//	 *								inputs are of different sizes or if less than
//	 *								three input samples are provided.
//	 */
//	/* * КОНСТРУКТОР *********************************************** ************** /
//
//
//        / **
//        * Выполняет преобразование Фурье и сохраняет реальные и мнимые результаты.
//* Входные сигналы дополняются нулями, если их длина не равна
//* сила 2.
//        *
//        * @param real_input Реальная часть сигнала, подлежащего преобразованию.
//* @param imaary_input Мнимая часть сигнала, которая будет.
//* Трансформируются. Это может быть нулевым, если сигнал
//* вполне реально.
//        * @param inverse_transform Значение false подразумевает, что форвард
//* должно применяться преобразование, а значение
//* true означает, что обратное преобразование является tob
//* быть примененным.
//        * @param use_hanning_window Значение true означает, что окно Hanning
//* будет применяться к real_input. Ценность
//* Вальс приведет к применению
//* окно Ханнинга.
//        * @throws Exception Выдает исключение, если реальное и воображаемое
//* Входные данные имеют разные размеры или если меньше
//* предоставляются три входных образца.
//        * /
//        * */
//
//
//	public FFT( double[] real_input,
//	            double[] imaginary_input,
//	            boolean inverse_transform,
//	            boolean use_hanning_window )
//		throws Exception
//	{
//
//		// Throw an exception if non-matching input signals are provided
//        //// Создаем исключение, если не совпадают входные сигналы
//		if (imaginary_input != null)
//			if (real_input.length != imaginary_input.length)
//				throw new Exception("Imaginary and real inputs are of different sizes.");
//
//		// Throw an exception if less than three samples are provided
//        // Создаем исключение, если предоставлено менее трех образцов
//		if (real_input.length < 3)
//			throw new Exception( "Only " + real_input.length + " samples provided.\n" +
//			                     "At least three are needed." );
//
//		// Verify that the input size has a number of samples that is a
//		// power of 2. If not, then increase the size of the array using
//		// zero-padding. Also creates a zero filled imaginary component
//		// of the input if none was specified.
//        //
//
//
//// Убедитесь, что входной размер имеет количество выборок, которое является
//// степень 2. Если нет, то увеличьте размер массива, используя
//// заполнение нулями. Также создает заполненный нулями мнимый компонент
//// ввода, если ничего не было указано.
//
//		int valid_size = ensureIsPowerOfN(real_input.length, 2);
//		if (valid_size != real_input.length)
//		{
//			double[] temp = new double[valid_size];
//			for (int i = 0; i < real_input.length; i++)
//				temp[i] = real_input[i];
//			for (int i = real_input.length; i < valid_size; i++)
//				temp[i] = 0.0;
//			real_input = temp;
//
//			if (imaginary_input == null)
//			{
//				imaginary_input = new double[valid_size];
//				for (int i = 0; i < imaginary_input.length; i++)
//					imaginary_input[i] = 0.0;
//			}
//			else
//			{
//				temp = new double[valid_size];
//				for (int i = 0; i < imaginary_input.length; i++)
//					temp[i] = imaginary_input[i];
//				for (int i = imaginary_input.length; i < valid_size; i++)
//					temp[i] = 0.0;
//				imaginary_input = temp;
//			}
//		}
//		else if (imaginary_input == null)
//		{
//			imaginary_input = new double[valid_size];
//			for (int i = 0; i < imaginary_input.length; i++)
//				imaginary_input[i] = 0.0;
//		}
//
//		// Instantiate the arrays to hold the output and copy the input
//		// to them, since the algorithm used here is self-processing
//        // Создание экземпляров массивов для хранения вывода и копирование ввода
//        // для них, так как алгоритм, используемый здесь, является самообработкой
//		real_output = new double[valid_size];
//		System.arraycopy(real_input, 0, real_output, 0, valid_size);
//		imaginary_output = new double[valid_size];
//		System.arraycopy(imaginary_input, 0, imaginary_output, 0, valid_size);
//
//		// Apply a Hanning window to the real values if this option is
//		// selected / Применить окно Ханнинга к реальным значениям, если эта опция
//        //// выбрано
//		if (use_hanning_window)
//		{
//			for (int i = 0; i < real_output.length; i++)
//			{
//				double hanning = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / valid_size);
//				real_output[i] *= hanning;
//			}
//		}
//
//		// Determine whether this is a forward or inverse transform
//        //// Определяем, является ли это прямое или обратное преобразование
//		int forward_transform = 1;
//		if (inverse_transform)
//			forward_transform = -1;
//
//		// Reorder the input data into reverse binary order
//        // // Переупорядочиваем входные данные в обратном двоичном порядке
//		double scale = 1.0;
//		int j = 0;
//		for (int i = 0; i < valid_size; ++i)
//		{
//			if (j >= i)
//			{
//				double tempr = real_output[j] * scale;
//				double tempi = imaginary_output[j] * scale;
//				real_output[j] = real_output[i] * scale;
//				imaginary_output[j] = imaginary_output[i] * scale;
//				real_output[i] = tempr;
//				imaginary_output[i] = tempi;
//			}
//			int m = valid_size / 2;
//			while (m >= 1 && j >= m)
//			{
//				j -= m;
//				m /= 2;
//			}
//			j += m;
//		}
//
//		// Perform the spectral recombination stage by stage
//        // Выполнение спектральной рекомбинации поэтапно
//		int stage = 0;
//		int max_spectra_for_stage; // максимый спектр
//		int step_size;  /// размер шага
//		for( max_spectra_for_stage = 1, step_size = 2 * max_spectra_for_stage;
//		     max_spectra_for_stage < valid_size;
//			 max_spectra_for_stage = step_size, step_size = 2 * max_spectra_for_stage)
//		{
//			double delta_angle = forward_transform * Math.PI / max_spectra_for_stage;
//
//			// Loop once for each individual spectra
//			for (int spectra_count = 0; spectra_count < max_spectra_for_stage; ++spectra_count)
//			{
//				double angle = spectra_count * delta_angle;
//				double real_correction = Math.cos(angle);
//				double imag_correction = Math.sin(angle);
//
//				int right = 0;
//				for (int left = spectra_count; left < valid_size; left += step_size)
//				{
//					right = left + max_spectra_for_stage;
//					double temp_real = real_correction * real_output[right] -
//					                   imag_correction * imaginary_output[right];
//					double temp_imag = real_correction * imaginary_output[right] +
//					                   imag_correction * real_output[right];
//					real_output[right] = real_output[left] - temp_real;
//					imaginary_output[right] = imaginary_output[left] - temp_imag;
//					real_output[left] += temp_real;
//					imaginary_output[left] += temp_imag;
//				}
//			}
//			max_spectra_for_stage = step_size;
//		}
//
//		// Set the angle and magnitude to null originally
//		output_angle = null;
//		output_power = null;
//		output_magnitude = null;
//	}
//
//    public FFT(int i) {
//    }
//
//
//    /* PUBLIC METHODS **********************************************************/
//
//
//	/**
//	 * Returns the magnitudes spectrum. It only makes sense to call
//	 * this method if this object was instantiated as a forward Fourier
//	 * transform.
//	 *
//	 * <p>Only the left side of the spectrum is returned, as the folded
//	 * portion of the spectrum is redundant for the purpose of the magnitude
//	 * spectrum. This means that the bins only go up to half of the
//	 * sampling rate.
//	 *
//	 * @return	The magnitude of each frequency bin.
//	 */
//	/*
//	/ **
//* Возвращает спектр величин. Имеет смысл использовать
//* этот метод, если этот объект был создан как прямой Фурье
//* преобразование.
//*
//* <p> Возвращается только левая часть спектра в сложенном виде.
//* часть спектра избыточна для целей величины
//* спектр. Это означает, что бункеры доходят только до половины
//* частота выборки.
//*
//
//* /
//	 */
//	public double[] getMagnitudeSpectrum()
//	{
//		// Only calculate the magnitudes if they have not yet been calculated
//        //// Рассчитать величины, только если они еще не были рассчитаны
//		if (output_magnitude == null)
//		{
//		    //количество развернутых бунов= мнимая_входная
//			int number_unfolded_bins = imaginary_output.length / 2;
//			//выходная магнитуда
//			output_magnitude = new double[number_unfolded_bins];
//			for(int i = 0; i < output_magnitude.length; i++)
//				output_magnitude[i] = ( Math.sqrt(real_output[i] * real_output[i] + imaginary_output[i] * imaginary_output[i]) ) / real_output.length;
//		}
//
//		// Return the magnitudes
//		return output_magnitude;
//	}
//
//
//	/**
//	 * Returns the power spectrum. It only makes sense to call
//	 * this method if this object was instantiated as a forward Fourier
//	 * transform.
//	 *
//	 * <p>Only the left side of the spectrum is returned, as the folded
//	 * portion of the spectrum is redundant for the purpose of the power
//	 * spectrum. This means that the bins only go up to half of the
//	 * sampling rate.
//	 *
//	 * @return	The magnitude of each frequency bin.
//	 */
//	/*
//     **
//     * Возвращает спектр мощности. Имеет смысл звонить
//     * этот метод, если этот объект был создан как прямой Фурье
//     * преобразование.
//     *
//     * <p> Возвращается только левая часть спектра в сложенном виде.
//     * часть спектра избыточна для целей питания
//     * спектр. Это означает, что бункеры доходят только до половины
//     * частота выборки.
//     *
//     * @return Величина каждого частотного бина.
//     * /
//	 */
//	public double[] getPowerSpectrum()
//	{
//		// Only calculate the powers if they have not yet been calculated
//        //// Рассчитать мощность, только если они еще не были рассчитаны
//		if (output_power == null)
//		{
//			int number_unfolded_bins = imaginary_output.length / 2;
//			output_power = new double[number_unfolded_bins];
//			for(int i = 0; i < output_power.length; i++)
//				output_power[i] = (real_output[i] * real_output[i] + imaginary_output[i] * imaginary_output[i]) / real_output.length;
//		}
//
//		// Return the power
//		return output_power;
//	}
//
//
//	/**
//	 * Returns the phase angle for each frequency bin. It only makes sense to
//	 * call this method if this object was instantiated as a forward Fourier
//	 * transform.
//	 *
//	 * <p>Only the left side of the spectrum is returned, as the folded
//	 * portion of the spectrum is redundant for the purpose of the phase
//	 * angles. This means that the bins only go up to half of the
//	 * sampling rate.
//	 *
//	 * @return	The phase angle for each frequency bin in degrees.
//	 */
//	/*
//	/ **
//* Возвращает фазовый угол для каждого частотного бина. Это имеет смысл только
//* вызовите этот метод, если этот объект был создан как прямой Фурье
//* преобразование.
//*
//* <p> Возвращается только левая часть спектра в сложенном виде.
//* часть спектра избыточна для целей фазы
//* углы. Это означает, что бункеры доходят только до половины
//* частота выборки.
//*
//* @return Фазовый угол для каждого частотного бина в градусах.
//* /
//	 */
//	public double[] getPhaseAngles()
//	{
//		// Only calculate the angles if they have not yet been calculated
//		if (output_angle == null)
//		{
//			int number_unfolded_bins = imaginary_output.length / 2;
//			output_angle = new double[number_unfolded_bins];
//			for(int i = 0; i < output_angle.length; i++)
//			{
//				if(imaginary_output[i] == 0.0 && real_output[i] == 0.0)
//					output_angle[i] = 0.0;
//				else
//					output_angle[i] = Math.atan(imaginary_output[i] / real_output[i]) * 180.0 / Math.PI;
//
//				if(real_output[i] < 0.0 && imaginary_output[i] == 0.0)
//					output_angle[i] = 180.0;
//				else if(real_output[i] < 0.0 && imaginary_output[i] == -0.0)
//					output_angle[i] = -180.0;
//				else if(real_output[i] < 0.0 && imaginary_output[i] > 0.0)
//					output_angle[i] += 180.0;
//				else if(real_output[i] < 0.0 && imaginary_output[i] < 0.0)
//					output_angle[i] += -180.0;
//			}
//		}
//
//		// Return the phase angles
//		return output_angle;
//	}
//
//
//	/**
//	 * Returns the frequency bin labels for each bin referred to by the
//	 * real values, imaginary values, magnitudes and phase angles as
//	 * determined by the given sampling rate.
//	 *
//	 * @param	sampling_rate	The sampling rate that was used to perform
//	 *							the FFT.
//	 * @return					The bin labels.
//	 */
//	/*
//            * Возвращает метки частотных бинов для каждого бина, указанного
//* реальные значения, мнимые значения, величины и фазовые углы как
//* определяется заданной частотой дискретизации.
//        *
//        * @param sampling_rate Частота дискретизации, которая использовалась для выполнения
//* БПФ.
//* @return Метки для мусора.
//*
//* */
//
//	public double[] getBinLabels(double sampling_rate)
//	{
//		int number_bins = real_output.length;
//		double bin_width = sampling_rate / (double) number_bins;
//		int number_unfolded_bins = imaginary_output.length / 2;
//		double[] labels = new double[number_unfolded_bins];
//		labels[0] = 0.0;
//		for (int bin = 1; bin < labels.length; bin++)
//			labels[bin] = bin * bin_width;
//		return labels;
//	}
//
//
//	/**
//	 * Returns the real values as calculated by the FFT.
//	 ** Возвращает реальные значения, рассчитанные БПФ.
//	 * @return	The real values.
//	 */
//	public double[] getRealValues()
//	{
//		return real_output;
//	}
//
//
//	/**
//	 * Returns the real values as calculated by the FFT.
//	 ** Возвращает реальные значения, рассчитанные БПФ.
//	 * @return	The real values.
//	 */
//	public double[] getImaginaryValues()
//	{
//		return imaginary_output;
//	}
//
//	/**
//	 * If the given x is a power of the given n, then x is returned.
//	 * If not, then the next value above the given x that is a power
//	 * of n is returned.
//	 *
//	 * Both x and n must be greater than zero.
//	 *
//	 * @param	x	The value to ensure is a power of n.
//	 * @param	n	The power to base x's validation on.
//	 */
//	/*/ **
//* Если данный x является степенью данного n, то возвращается x.
//* Если нет, то следующее значение выше заданного x, которое является степенью
//* из n возвращается.
//*
//* И x, и n должны быть больше нуля.
//*
//* @param x Значение для обеспечения является степенью n.
//* @param n Мощность, на которой основывается проверка x.
//* */
//
//	private static int ensureIsPowerOfN(int x, int n)
//	{
//		double log_value = logBaseN((double) x, (double) n);
//		int log_int = (int) log_value;
//		int valid_size = pow(n, log_int);
//		if (valid_size != x)
//			valid_size = pow(n, log_int + 1);
//		return valid_size;
//	}
//
//
//	/**
//	 * Returns the logarithm of the specified base of the given number.
//	 *  Both x and n must be greater than zero.
//	 * @param	x	The value to find the log of.
//	 * @param	n	The base of the logarithm.
//     *              Возвращает логарифм указанной базы данного числа.
//     * *
//     * * И x, и n должны быть больше нуля.
//     * *
//     * * @param x Значение для поиска журнала.
//     * * @param n Основание логарифма.
//	 */
//	private static double logBaseN(double x, double n)
//	{
//		return (Math.log10(x) / Math.log10(n));
//	}
//
//
//	/**
//	 * Returns the given a raised to the power of the given b.
//	 *
//	 *  b must be greater than zero.
//	 *
//	 * @param	a	The base.
//	 * @param	b	The exponent.
//     *              * Возвращает данное повышение в силу данного б.
//     * *
//     * * b должно быть больше нуля.
//     * *
//     * * @param a База.
//     * * @param b Показатель степени.
//	 */
//	private static int pow(int a, int b)
//	{
//		int result = a;
//		for (int i = 1; i < b; i++)
//			result *= a;
//		return result;
//	}
//    public int getN() {
//        return n;
//    }
//
//    public int getM() {
//        return m;
//    }
//
//    public int getK() {
//        return k;
//    }
//}

private int n;
private int m;
private int k;
private float[] cos;
private float[] sin;

public FFT(int n) {
        this.m = (n == 0)? 0 : (32 - Integer.numberOfLeadingZeros(n - 1));
        this.n = n = 1 << m;
        this.k = n / 2 + 1;

        cos = new float[n / 2];
        sin = new float[n / 2];
        for (int i = 0; i < n / 2; i++) {
        cos[i] = (float) (Math.cos(-2 * Math.PI * i / n));
        sin[i] = (float) (Math.sin(-2 * Math.PI * i / n));
        }
        }

public int getN() {
        return n;
        }

public int getM() {
        return m;
        }

public int getK() {
        return k;
        }

private void reverse(float[] x, float[] y) {
        // TODO: bit reversal can only be disabled when compacting properly (please do this!)
        int i, j, n1, n2;
        float t1;

        j = 0;
        n2 = n / 2;
        for (i = 1; i < n - 1; i++) {
        n1 = n2;
        while (j >= n1) {
        j = j - n1;
        n1 = n1 / 2;
        }
        j = j + n1;

        if (i < j) {
        t1 = x[i];
        x[i] = x[j];
        x[j] = t1;
        t1 = y[i];
        y[i] = y[j];
        y[j] = t1;
        }
        }
        }

public void fft(float[] x, float[] y, boolean reverse) {
        if ((x.length != n) || (y.length != n)) {
        throw new RuntimeException();
        }

        int i, j, k, n1, n2, a;
        float c, s, t1, t2;

        if(reverse) {
        reverse(x, y);
        }

        // FFT
        n1 = 0;
        n2 = 1;

        for (i = 0; i < m; i++) {
        n1 = n2;
        n2 = n2 + n2;
        a = 0;

        for (j = 0; j < n1; j++) {
        c = cos[a];
        s = sin[a];
        a += 1 << (m - i - 1);

        for (k = j; k < n; k = k + n2) {
        t1 = c * x[k + n1] - s * y[k + n1];
        t2 = s * x[k + n1] + c * y[k + n1];
        x[k + n1] = x[k] - t1;
        y[k + n1] = y[k] - t2;
        x[k] = x[k] + t1;
        y[k] = y[k] + t2;
        }
        }
        }
        }

public void ifft(float[] x, float[] y, boolean reverse) {
        if ((x.length != n) || (y.length != n)) {
        throw new RuntimeException();
        }

        int i, j, k, n1, n2, a;
        float c, s, t1, t2;

        if(reverse) {
        reverse(x, y);
        }

        // FFT
        n1 = 0;
        n2 = 1;

        for (i = 0; i < m; i++) {
        n1 = n2;
        n2 = n2 + n2;
        a = 0;

        for (j = 0; j < n1; j++) {
        c = cos[a];
        s = sin[a];
        a += 1 << (m - i - 1);

        for (k = j; k < n; k = k + n2) {
        t1 = c * x[k + n1] + s * y[k + n1];
        t2 = s * x[k + n1] - c * y[k + n1];
        x[k + n1] = x[k] - t1;
        y[k + n1] = y[k] + t2;
        x[k] = x[k] + t1;
        y[k] = y[k] - t2;
        }
        }
        }
        }

public void scale(float[] x, float[] y) {
        if ((x.length != n) || (y.length != n)) {
        throw new RuntimeException();
        }

        for (int i = 0; i < n; i++) {
        x[i] /= n;
        y[i] /= n;
        }
        }
        }