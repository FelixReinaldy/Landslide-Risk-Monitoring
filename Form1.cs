using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using System.Drawing;

namespace LandslideFusion
{
    public partial class Form1 : Form
    {
        // ================== SAMPLING GLOBAL ==================
        // Frekuensi sampling dikunci di 100 Hz untuk semua sensor
        private const double Fs = 100.0;         // Hz
        private const double Dt = 1.0 / Fs;      // 0.01 s

        // Model orde-1 gabungan: τ dR/dt + R = u(t)
        // τ bisa kamu jelaskan sebagai parameter desain / dari literatur
        private const double TauModel = 5.0;     // detik

        // buffer sinyal per-sensor
        private List<double> soilMoist = new();
        private List<double> slopeDeg = new();
        private List<double> accel = new();     // ground acceleration geophone [m/s^2]
        private List<double> pwp = new();
        private List<double> rainInt = new();

        private int N = 2048;
        private double dtSeconds = Dt;          // pastikan selalu 0.01 s
        private bool _ready = false;

        public Form1()
        {
            InitializeComponent();
            AutoScaleMode = AutoScaleMode.Font;
            panelCharts.AutoScroll = true;
            Load += Form1_Load;
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            SafeInitCharts();
            TryLoadBackground();

            _ready = true;
            RecomputeAll();
        }

        private void TryLoadBackground()
        {
            try
            {
                string fileName = "landslide_bg.png";
                string imgPath = System.IO.Path.Combine(Application.StartupPath, fileName);

                if (System.IO.File.Exists(imgPath))
                {
                    panel3D.BackgroundImage = Image.FromFile(imgPath);
                    panel3D.BackgroundImageLayout = ImageLayout.Stretch;
                }
            }
            catch { }
        }

        /* ================== EVENT INPUT ================== */
        private void AnyInputChanged(object sender, EventArgs e)
        {
            if (!_ready) return;

            // dtSeconds TIDAK lagi bisa diubah dari UI, selalu Dt = 0.01 s (Fs = 100 Hz)

            lblRAW.Text = $"RAW: {tbRAW.Value}";
            lblVdiff.Text = $"Vout-V0: {tbVdiff.Value / 100.0:0.00} V,  k={numK.Value:0.###} V/deg";
            lblVacc.Text = $"V: {tbVacc.Value / 100.0:0.00} V";
            lblF.Text = $"f: {tbF.Value} Hz,  f₀: {tbF0.Value} Hz,  a={numKp.Value:E2}";
            lblDm.Text = $"Δm/step: {tbDm.Value / 10.0:0.0} g,  A={numArea.Value:0.###} m²";

            RecomputeAll();
        }

        private void btnRumusMoist_Click(object _, EventArgs __) =>
            txtFormula.Text = "TEROS-12: VWC = 3.879 × 10⁻⁴ × RAW − 0.6956";

        private void btnRumusSlope_Click(object _, EventArgs __) =>
            txtFormula.Text = "A801: θ = (Vout − V₀) / k";

        // ===> RUMUS GEOPHONE 4.5 Hz
        private void btnRumusVib_Click(object _, EventArgs __) =>
            txtFormula.Text = "Geophone 4.5 Hz: v_g = V / 82  [m/s], a_g ≈ dv_g/dt";

        // ===> RUMUS PORE PRESSURE BERDASARKAN JURNAL
        private void btnRumusPress_Click(object _, EventArgs __) =>
            txtFormula.Text = "VWP-3000: P = a × (f² − f₀²)";

        private void btnRumusRain_Click(object _, EventArgs __) =>
            txtFormula.Text = "Pluvio²: R = 3.6 × Δm / (A × Δt) [g,m²,detik]";

        /* ================== PERHITUNGAN DAN PLOT ================== */
        private void RecomputeAll()
        {
            // dtSeconds adalah Dt (= 0.01 s, fs = 100 Hz)
            double dt = dtSeconds;

            // --- simulasi sinyal raw per sensor ---
            var rawSig = Const(N, (double)tbRAW.Value);
            rawSig = rawSig.Select((v, i) => v + 30 * Math.Sin(2 * Math.PI * i * dt / 200.0)).ToList();

            var vdiffSig = Sin(N, dt, 1.0 / 150.0, tbVdiff.Value / 100.0, 0.0);

            // Tegangan geophone (V), 0–2 V dari trackbar
            var vSig = Sin(N, dt, 1.0 / 10.0, tbVacc.Value / 100.0, 0.0);

            var fSig = Const(N, (double)tbF.Value)
                       .Select((v, i) => v + 5 * Math.Sin(2 * Math.PI * i * dt / 300.0)).ToList();
            double f0 = tbF0.Value;

            var dmStep = Const(N, tbDm.Value / 10.0);

            // --- konversi ke besaran fisik ---

            // TEROS-12: VWC
            soilMoist = rawSig.Select(RAW_to_VWC).ToList();

            // A801 — θ = (Vout − V0) / k
            double k = (double)numK.Value;
            slopeDeg = vdiffSig.Select(x => x / k).ToList();

            // ================= GEOPHONE 4.5 Hz =================
            // Datasheet: sensitivitas S ≈ 82 V/(m/s)
            // Langkah 1: v_g = V / S              → ground velocity [m/s]
            // Langkah 2: a_g = dv_g/dt (diskrit)  → ground acceleration [m/s^2]
            const double S_geophone = 82.0;   // V per (m/s)

            // konversi tegangan geophone menjadi ground velocity
            var vel = vSig.Select(v => v / S_geophone).ToList();   // [m/s]

            // diferensial numerik: velocity -> acceleration
            accel = VelocityToAcceleration(vel, dt);               // [m/s^2]

            // =============== VWP-3000 PORE PRESSURE ===============
            // Model berdasar jurnal: P = a (f^2 - f0^2)
            double aCoeff = (double)numKp.Value;   // koefisien kalibrasi a [pressure per Hz^2]
            pwp = fSig.Select(x => aCoeff * (x * x - f0 * f0)).ToList();

            // Pluvio
            rainInt = RainIntensityFromDm(dmStep, (double)numArea.Value, dt);

            // --- Plot time & frequency per sensor ---
            PlotTimeAndFreq(chTimeMoist, chFreqMoist, soilMoist, dt, "VWC (m³/m³)");
            PlotTimeAndFreq(chTimeSlope, chFreqSlope, slopeDeg, dt, "Tilt (deg)");
            // sekarang menampilkan ground acceleration geophone
            PlotTimeAndFreq(chTimeAcc, chFreqAcc, accel, dt, "Ground Acceleration (m/s²)");
            PlotTimeAndFreq(chTimePwp, chFreqPwp, pwp, dt, "Pore Pressure");
            PlotTimeAndFreq(chTimeRain, chFreqRain, rainInt, dt, "Rain (mm/h)");

            // ====== GABUNGAN 5 SENSOR: DIPAKAI HANYA UNTUK ODE, S-DOMAIN, Z-DOMAIN ======
            // u(t) = w1 z1 + ... + w5 z5  (gabungan statis dari sensor)
            var u = BuildFusionInput();

            // R(t) = solusi diskrit ODE orde-1: τ dR/dt + R = u(t)
            var R = SolveFirstOrderResponse(u, dtSeconds, TauModel);

            // Teks di UI: definisi z_i, u(t), dan persamaan diferensial gabungan
            txtRumusGabungan.Text =
                "z₁ = z(3.879×10⁻⁴·RAW − 0.6956)" + Environment.NewLine +
                "z₂ = z((Vout − V₀)/k)" + Environment.NewLine +
                "z₃ = z(a_g)" + Environment.NewLine +
                "z₄ = z(a·(f² − f₀²))" + Environment.NewLine +
                "z₅ = z(3.6·Δm/(A·Δt))" + Environment.NewLine +
                "u(t) = w₁·z₁ + w₂·z₂ + w₃·z₃ + w₄·z₄ + w₅·z₅" + Environment.NewLine +
                "Persamaan diferensial gabungan:" + Environment.NewLine +
                $"τ·dR/dt + R(t) = u(t),  τ = {TauModel:0.###} s";

            // Plot S-domain dan Z-domain berdasarkan model orde-1 dengan τ tetap
            PlotSAndZ(TauModel);
        }

        // u(t): gabungan statis dari 5 sensor (weighted z-score)
        private List<double> BuildFusionInput()
        {
            int n = new[] {
                soilMoist.Count, slopeDeg.Count, accel.Count, pwp.Count, rainInt.Count }.Min();

            const double wMoist = 0.25;
            const double wSlope = 0.25;
            const double wAccel = 0.20;
            const double wPwp = 0.15;
            const double wRain = 0.15;

            var u = new List<double>(n);
            for (int i = 0; i < n; i++)
            {
                double s = 0;
                s += wMoist * Z(soilMoist, i, n);
                s += wSlope * Z(slopeDeg, i, n);
                s += wAccel * Z(accel, i, n);
                s += wPwp * Z(pwp, i, n);
                s += wRain * Z(rainInt, i, n);
                u.Add(s);
            }
            return u;
        }

        // R(t): solusi diskrit dari ODE orde-1: τ dR/dt + R = u(t)
        // Diskret: R[k] = a·R[k-1] + b·u[k], dengan a = e^{-T/τ}, b = 1 - a
        private List<double> SolveFirstOrderResponse(List<double> u, double dt, double tau)
        {
            var R = new List<double>(u.Count);
            if (u.Count == 0) return R;

            double a = Math.Exp(-dt / tau);
            double b = 1.0 - a;

            double Rk = u[0];  // kondisi awal: R(0) ≈ u(0)
            R.Add(Rk);

            for (int k = 1; k < u.Count; k++)
            {
                Rk = a * Rk + b * u[k];
                R.Add(Rk);
            }

            return R;
        }

        private void PlotSAndZ(double tau)
        {
            // ===== S-DOMAIN =====
            // G(s) = 1 / (τ s + 1)
            txtSDomain.Text = $"G(s) = 1 / ({tau:0.###}·s + 1)";

            SetupPoleZeroChart(chartSDomain, false);
            double sp = -1.0 / tau; // pole di s-plane
            chartSDomain.Series["Poles"].Points.AddXY(sp, 0.0);

            // ===== Z-DOMAIN (ZOH) =====
            // G(z) = (1 - e^{-T/τ}) z^{-1} / (1 - e^{-T/τ} z^{-1})
            // Pole: z_p = e^{-T/τ}, Zero: z = 0
            SetupPoleZeroChart(chartZDomain, true);

            double T = dtSeconds;
            double zp = Math.Exp(-T / tau);  // real, di dalam unit circle
            double zz = 0.0;                 // zero di origin

            chartZDomain.Series["Poles"].Points.AddXY(zp, 0.0);
            chartZDomain.Series["Zeros"].Points.AddXY(zz, 0.0);
        }

        /* ================== PLOT HELPERS ================== */

        private void PlotTimeAndFreq(Chart chTime, Chart chFreq, List<double> x, double dt, string yLabel)
        {
            EnsureChartSetup(chTime);
            EnsureChartSetup(chFreq);

            if (x == null || x.Count == 0) return;

            // --- time domain ---
            chTime.Series[0].Points.Clear();
            for (int i = 0; i < x.Count; i++)
                chTime.Series[0].Points.AddXY(i * dt, x[i]);

            chTime.ChartAreas[0].AxisX.Title = "Time (s)";
            chTime.ChartAreas[0].AxisY.Title = yLabel;

            // --- frequency domain ---
            var spec = FFTMag(x);
            if (spec.Count < 2) return;

            double fs = 1.0 / dt;
            int half = spec.Count / 2;
            chFreq.Series[0].Points.Clear();

            for (int k = 1; k < half; k++)
            {
                double f = k * fs / spec.Count;
                chFreq.Series[0].Points.AddXY(f, spec[k]);
            }

            chFreq.ChartAreas[0].AxisX.Title = "Frequency (Hz)";
            chFreq.ChartAreas[0].AxisY.Title = "Magnitude";
        }

        private void SafeInitCharts()
        {
            foreach (var ch in GetAllCharts(this))
                EnsureChartSetup(ch);

            SetupPoleZeroChart(chartSDomain, false);
            SetupPoleZeroChart(chartZDomain, true);
        }

        private IEnumerable<Chart> GetAllCharts(Control root)
        {
            foreach (Control c in root.Controls)
            {
                if (c is Chart ch)
                    yield return ch;

                foreach (var sub in GetAllCharts(c))
                    yield return sub;
            }
        }

        private void EnsureChartSetup(Chart c)
        {
            if (c.ChartAreas.Count == 0)
            {
                var ca = new ChartArea("ca");
                ca.AxisX.MajorGrid.Enabled = false;
                ca.AxisY.MajorGrid.Enabled = true;
                c.ChartAreas.Add(ca);
            }

            if (c.Series.Count == 0)
            {
                var s = new Series("data")
                { ChartType = SeriesChartType.Line, BorderWidth = 2 };
                c.Series.Add(s);
            }
        }

        private void SetupPoleZeroChart(Chart chart, bool zPlane)
        {
            chart.ChartAreas.Clear();
            var ca = new ChartArea("pz");

            ca.AxisX.Crossing = 0;
            ca.AxisY.Crossing = 0;
            ca.AxisX.MajorGrid.Enabled = false;
            ca.AxisY.MajorGrid.Enabled = false;

            chart.ChartAreas.Add(ca);
            chart.Series.Clear();

            if (zPlane)
            {
                ca.AxisX.Minimum = -1.5;
                ca.AxisX.Maximum = 1.5;
                ca.AxisY.Minimum = -1.5;
                ca.AxisY.Maximum = 1.5;
            }
            else
            {
                ca.AxisX.Minimum = -5;
                ca.AxisX.Maximum = 1;
                ca.AxisY.Minimum = -5;
                ca.AxisY.Maximum = 5;

                var shading = new Series("UnstableRegion")
                {
                    ChartType = SeriesChartType.Area,
                    Color = Color.FromArgb(60, Color.Red),
                    BorderWidth = 0
                };

                shading.Points.AddXY(0.0, ca.AxisY.Minimum);
                shading.Points.AddXY(ca.AxisX.Maximum, ca.AxisY.Minimum);
                shading.Points.AddXY(ca.AxisX.Maximum, ca.AxisY.Maximum);
                shading.Points.AddXY(0.0, ca.AxisY.Maximum);

                chart.Series.Add(shading);
            }

            var pole = new Series("Poles")
            {
                ChartType = SeriesChartType.Point,
                MarkerStyle = MarkerStyle.Cross,
                MarkerSize = 12
            };

            var zero = new Series("Zeros")
            {
                ChartType = SeriesChartType.Point,
                MarkerStyle = MarkerStyle.Circle,
                MarkerSize = 8
            };

            chart.Series.Add(pole);
            chart.Series.Add(zero);

            if (zPlane)
            {
                var uc = new Series("UnitCircle")
                {
                    ChartType = SeriesChartType.Line,
                    BorderDashStyle = ChartDashStyle.Dash
                };

                for (int i = 0; i <= 256; i++)
                {
                    double a = 2 * Math.PI * i / 256.0;
                    uc.Points.AddXY(Math.Cos(a), Math.Sin(a));
                }

                chart.Series.Add(uc);
            }
        }

        /* ================== UTILITIES & RUMUS ================== */

        private static double RAW_to_VWC(double raw)
            => 3.879e-4 * raw - 0.6956;

        private static List<double> RainIntensityFromDm(List<double> dm_g_step, double area_m2, double dt_s)
        {
            var cum = Cumulative(dm_g_step);
            var rate_g_per_s = new List<double> { 0 };

            for (int i = 1; i < cum.Count; i++)
                rate_g_per_s.Add((cum[i] - cum[i - 1]) / dt_s);

            return rate_g_per_s.Select(r => 3.6 * r / area_m2).ToList();
        }

        // konversi velocity [m/s] menjadi acceleration [m/s^2] dengan diferensial diskrit
        private static List<double> VelocityToAcceleration(List<double> vel, double dt)
        {
            var acc = new List<double>(vel.Count);

            if (vel.Count == 0)
                return acc;

            // asumsi awal: percepatan awal = 0
            acc.Add(0.0);

            for (int i = 1; i < vel.Count; i++)
            {
                double a = (vel[i] - vel[i - 1]) / dt;   // m/s^2
                acc.Add(a);
            }

            return acc;
        }

        private static IEnumerable<double> Sin(int N, double dt, double f, double amp, double noiseAmp)
        {
            var rnd = new Random(0);
            for (int i = 0; i < N; i++)
            {
                double n = noiseAmp * (rnd.NextDouble() - 0.5);
                yield return amp * Math.Sin(2 * Math.PI * f * i * dt) + n;
            }
        }

        private static List<double> Const(int N, double v)
        {
            var y = new List<double>(N);
            for (int i = 0; i < N; i++)
                y.Add(v);
            return y;
        }

        private static List<double> Cumulative(List<double> step)
        {
            double a = 0;
            var y = new List<double>(step.Count);
            foreach (var s in step)
            {
                a += s;
                y.Add(a);
            }
            return y;
        }

        private static double Z(List<double> x, int idx, int N)
        {
            double m = x.Take(N).Average();
            double s = Math.Sqrt(x.Take(N)
                         .Select(v => (v - m) * (v - m)).Average() + 1e-12);

            return (x[idx] - m) / s;
        }

        private static List<double> FFTMag(List<double> x)
        {
            int n = 1;
            while (n < x.Count) n <<= 1;

            Complex[] a = new Complex[n];

            for (int i = 0; i < x.Count; i++)
                a[i] = new Complex(x[i], 0);

            for (int i = x.Count; i < n; i++)
                a[i] = Complex.Zero;

            FFT(a);

            return a.Select(c => c.Magnitude).ToList();
        }

        private static void FFT(Complex[] a)
        {
            int n = a.Length;

            for (int i = 1, j = 0; i < n; i++)
            {
                int bit = n >> 1;
                for (; (j & bit) != 0; bit >>= 1)
                    j ^= bit;
                j ^= bit;
                if (i < j)
                    (a[i], a[j]) = (a[j], a[i]);
            }

            for (int len = 2; len <= n; len <<= 1)
            {
                double ang = 2 * Math.PI / len;
                Complex wlen = new Complex(Math.Cos(ang), Math.Sin(ang));

                for (int i = 0; i < n; i += len)
                {
                    Complex w = Complex.One;

                    for (int j = 0; j < len / 2; j++)
                    {
                        Complex u = a[i + j];
                        Complex v = a[i + j + len / 2] * w;

                        a[i + j] = u + v;
                        a[i + j + len / 2] = u - v;

                        w *= wlen;
                    }
                }
            }
        }
    }
}
