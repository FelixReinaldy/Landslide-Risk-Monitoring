using System.Windows.Forms;

namespace LandslideFusion
{
    partial class Form1
    {
        private System.ComponentModel.IContainer components = null;
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null)) components.Dispose();
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code
        private void InitializeComponent()
        {
            // === PANEL KIRI (input sensor) ===
            panelLeft = new Panel { Dock = DockStyle.Left, Width = 300 };

            // 1) TEROS-12
            groupMoist = new GroupBox { Text = "Soil Moisture (TEROS-12)", Dock = DockStyle.Top, Height = 120 };
            btnRumusMoist = new Button { Text = "Rumus", Left = 15, Top = 25, Width = 90 };
            btnRumusMoist.Click += btnRumusMoist_Click;

            tbRAW = new TrackBar
            {
                Left = 120,
                Top = 20,
                Width = 150,
                Minimum = 0,
                Maximum = 4095,
                Value = 2000,
                TickFrequency = 256
            };
            tbRAW.ValueChanged += AnyInputChanged;

            lblRAW = new Label
            {
                Left = 120,
                Top = 80,
                Width = 160,
                Text = "RAW: 2000"
            };

            groupMoist.Controls.AddRange(new Control[] { btnRumusMoist, tbRAW, lblRAW });

            // 2) A801 – pakai k (V/deg)
            groupSlope = new GroupBox
            {
                Text = "Slope (A801)",
                Dock = DockStyle.Top,
                Height = 140
            };

            btnRumusSlope = new Button
            {
                Text = "Rumus",
                Left = 15,
                Top = 25,
                Width = 90
            };
            btnRumusSlope.Click += btnRumusSlope_Click;

            tbVdiff = new TrackBar
            {
                Left = 120,
                Top = 20,
                Width = 150,
                Minimum = -200,
                Maximum = 200,
                Value = 100,
                TickFrequency = 50
            };
            tbVdiff.ValueChanged += AnyInputChanged;

            // === NUMERIC K ===
            numK = new NumericUpDown
            {
                Left = 120,
                Top = 70,
                Width = 80,
                DecimalPlaces = 4,
                Increment = 0.001M,
                Minimum = 0.001M,
                Maximum = 10.000M,
                Value = 0.050M        // default 50 mV/deg
            };
            numK.ValueChanged += AnyInputChanged;

            lblVdiff = new Label
            {
                Left = 120,
                Top = 100,
                Width = 160,
                Text = "Vout-V0: 1.00 V,  k=0.050 V/deg"
            };

            groupSlope.Controls.AddRange(new Control[] {
                btnRumusSlope, tbVdiff, numK, lblVdiff
            });

            // 3) Geophone 4.5 Hz
            groupVib = new GroupBox
            {
                Text = "Vibration (Geophone 4.5 Hz)",
                Dock = DockStyle.Top,
                Height = 120
            };

            btnRumusVib = new Button { Text = "Rumus", Left = 15, Top = 25, Width = 90 };
            btnRumusVib.Click += btnRumusVib_Click;

            tbVacc = new TrackBar
            {
                Left = 120,
                Top = 20,
                Width = 150,
                Minimum = 0,
                Maximum = 200,
                Value = 100,
                TickFrequency = 25
            };
            tbVacc.ValueChanged += AnyInputChanged;

            lblVacc = new Label
            {
                Left = 120,
                Top = 80,
                Width = 160,
                Text = "V: 1.00 V"
            };

            groupVib.Controls.AddRange(new Control[] {
                btnRumusVib, tbVacc, lblVacc
            });

            // 4) VWP-3000
            groupPress = new GroupBox
            {
                Text = "Pore Pressure (VWP-3000)",
                Dock = DockStyle.Top,
                Height = 230
            };
            btnRumusPress = new Button
            {
                Text = "Rumus",
                Left = 15,
                Top = 25,
                Width = 90
            };
            btnRumusPress.Click += btnRumusPress_Click;

            numKp = new NumericUpDown
            {
                Left = 120,
                Top = 25,
                Width = 120,
                DecimalPlaces = 6,
                Increment = 0.000001M,
                Minimum = -1000,
                Maximum = 1000,
                Value = 0.000001M
            };
            numKp.ValueChanged += AnyInputChanged;

            var lblKp = new Label
            {
                Left = 250,
                Top = 28,
                AutoSize = true,
                Text = "a"
            };

            var lblFcap = new Label
            {
                Left = 120,
                Top = 58,
                AutoSize = true,
                Text = "f (Hz)"
            };

            tbF = new TrackBar
            {
                Left = 120,
                Top = 75,
                Width = 200,
                Minimum = 0,
                Maximum = 2000,
                Value = 1000,
                TickStyle = TickStyle.None
            };
            tbF.ValueChanged += AnyInputChanged;

            var lblF0cap = new Label
            {
                Left = 120,
                Top = 105,
                AutoSize = true,
                Text = "f₀ (Hz)"
            };

            tbF0 = new TrackBar
            {
                Left = 120,
                Top = 122,
                Width = 200,
                Minimum = 0,
                Maximum = 2000,
                Value = 900,
                TickStyle = TickStyle.None
            };
            tbF0.ValueChanged += AnyInputChanged;

            lblF = new Label
            {
                Left = 15,
                Top = 180,
                AutoSize = true,
                Text = "f: 1000 Hz, f₀: 900 Hz, a=1.0E-06"
            };

            groupPress.Controls.AddRange(new Control[] {
                btnRumusPress, numKp, lblKp,
                lblFcap, tbF, lblF0cap, tbF0, lblF
            });

            // 5) Pluvio²
            groupRain = new GroupBox
            {
                Text = "Rain Intensity (Pluvio²)",
                Dock = DockStyle.Top,
                Height = 150
            };
            btnRumusRain = new Button { Text = "Rumus", Left = 15, Top = 25, Width = 90 };
            btnRumusRain.Click += btnRumusRain_Click;

            tbDm = new TrackBar
            {
                Left = 120,
                Top = 20,
                Width = 150,
                Minimum = 0,
                Maximum = 500,
                Value = 10,
                TickFrequency = 50
            };
            tbDm.ValueChanged += AnyInputChanged;

            numArea = new NumericUpDown
            {
                Left = 120,
                Top = 70,
                Width = 80,
                DecimalPlaces = 3,
                Increment = 0.001M,
                Minimum = 0.001M,
                Maximum = 2.000M,
                Value = 0.020M
            };
            numArea.ValueChanged += AnyInputChanged;

            lblDm = new Label
            {
                Left = 120,
                Top = 100,
                Width = 170,
                Text = "Δm/step: 1.0 g, A=0.020 m²"
            };

            groupRain.Controls.AddRange(new Control[] {
                btnRumusRain, tbDm, numArea, lblDm
            });

            panelLeft.Controls.AddRange(new Control[] {
                groupRain, groupPress, groupVib, groupSlope, groupMoist
            });

            // === PANEL KANAN (info) ===
            panelRight = new Panel { Dock = DockStyle.Right, Width = 300 };

            label1 = new Label
            {
                Text = "Info / Rumus",
                Left = 20,
                Top = 20,
                AutoSize = true
            };

            txtFormula = new TextBox
            {
                Left = 20,
                Top = 50,
                Width = 260,
                Height = 160,
                Multiline = true,
                ReadOnly = true
            };

            panelRight.Controls.AddRange(new Control[] {
                label1, txtFormula
            });

            // === PANEL TENGAH ===
            panelCenter = new Panel { Dock = DockStyle.Fill };
            panelCharts = new Panel { Dock = DockStyle.Fill, AutoScroll = true };

            tlp = new TableLayoutPanel
            {
                Dock = DockStyle.Fill,
                ColumnCount = 2,
                RowCount = 5
            };
            tlp.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 50));
            tlp.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 50));
            for (int r = 0; r < 5; r++)
                tlp.RowStyles.Add(new RowStyle(SizeType.Percent, 20));

            chTimeMoist = new System.Windows.Forms.DataVisualization.Charting.Chart { Dock = DockStyle.Fill };
            chFreqMoist = new System.Windows.Forms.DataVisualization.Charting.Chart { Dock = DockStyle.Fill };

            chTimeSlope = new System.Windows.Forms.DataVisualization.Charting.Chart { Dock = DockStyle.Fill };
            chFreqSlope = new System.Windows.Forms.DataVisualization.Charting.Chart { Dock = DockStyle.Fill };

            chTimeAcc = new System.Windows.Forms.DataVisualization.Charting.Chart { Dock = DockStyle.Fill };
            chFreqAcc = new System.Windows.Forms.DataVisualization.Charting.Chart { Dock = DockStyle.Fill };

            chTimePwp = new System.Windows.Forms.DataVisualization.Charting.Chart { Dock = DockStyle.Fill };
            chFreqPwp = new System.Windows.Forms.DataVisualization.Charting.Chart { Dock = DockStyle.Fill };

            chTimeRain = new System.Windows.Forms.DataVisualization.Charting.Chart { Dock = DockStyle.Fill };
            chFreqRain = new System.Windows.Forms.DataVisualization.Charting.Chart { Dock = DockStyle.Fill };

            tlp.Controls.Add(chTimeMoist, 0, 0);
            tlp.Controls.Add(chFreqMoist, 1, 0);

            tlp.Controls.Add(chTimeSlope, 0, 1);
            tlp.Controls.Add(chFreqSlope, 1, 1);

            tlp.Controls.Add(chTimeAcc, 0, 2);
            tlp.Controls.Add(chFreqAcc, 1, 2);

            tlp.Controls.Add(chTimePwp, 0, 3);
            tlp.Controls.Add(chFreqPwp, 1, 3);

            tlp.Controls.Add(chTimeRain, 0, 4);
            tlp.Controls.Add(chFreqRain, 1, 4);

            panelCharts.Controls.Add(tlp);

            // === PANEL KANAN KECIL (3D + rumus gabungan + S/Z) ===
            panelRightCol = new Panel { Dock = DockStyle.Right, Width = 380 };

            lbl3D = new Label
            {
                Text = "Gambar 3D (placeholder)",
                Left = 10,
                Top = 0,
                AutoSize = true
            };

            panel3D = new Panel
            {
                Left = 10,
                Top = 20,
                Width = 360,
                Height = 160,
                BorderStyle = BorderStyle.FixedSingle
            };

            label4 = new Label
            {
                Text = "Rumus Gabungan",
                Left = 10,
                Top = 190,
                AutoSize = true
            };

            txtRumusGabungan = new TextBox
            {
                Left = 10,
                Top = 210,
                Width = 360,
                Height = 60,
                Multiline = true,
                ReadOnly = true
            };

            label5 = new Label
            {
                Text = "S Domain (pole)  |  Z Domain (pole/zero)",
                Left = 10,
                Top = 275,
                AutoSize = true
            };

            txtSDomain = new TextBox
            {
                Left = 10,
                Top = 295,
                Width = 360
            };

            chartSDomain = new System.Windows.Forms.DataVisualization.Charting.Chart
            {
                Left = 10,
                Top = 325,
                Width = 360,
                Height = 180
            };

            chartZDomain = new System.Windows.Forms.DataVisualization.Charting.Chart
            {
                Left = 10,
                Top = 510,
                Width = 360,
                Height = 200
            };

            panelRightCol.Controls.AddRange(new Control[] {
                lbl3D, panel3D,
                label4, txtRumusGabungan,
                label5, txtSDomain,
                chartSDomain, chartZDomain
            });

            panelCenter.Controls.Add(panelCharts);
            panelCenter.Controls.Add(panelRightCol);

            // === FORM ===
            Text = "Landslide Risk Monitoring - Sensor Fusion";
            ClientSize = new System.Drawing.Size(1360, 760);

            Controls.Add(panelCenter);
            Controls.Add(panelRight);
            Controls.Add(panelLeft);
        }
        #endregion

        // ==== FIELDS ====
        private Panel panelLeft, panelRight, panelCenter, panelCharts, panelRightCol, panel3D;

        private GroupBox groupMoist, groupSlope, groupVib, groupPress, groupRain;
        private Button btnRumusMoist, btnRumusSlope, btnRumusVib, btnRumusPress, btnRumusRain;

        private TrackBar tbRAW, tbVdiff, tbVacc, tbF, tbF0, tbDm;

        private NumericUpDown numK, numKp, numArea;

        private Label lblRAW, lblVdiff, lblVacc, lblF, lblDm, label1, lbl3D, label4, label5;

        private TextBox txtFormula, txtRumusGabungan, txtSDomain;

        private TableLayoutPanel tlp;

        private System.Windows.Forms.DataVisualization.Charting.Chart
            chTimeMoist, chFreqMoist,
            chTimeSlope, chFreqSlope,
            chTimeAcc, chFreqAcc,
            chTimePwp, chFreqPwp,
            chTimeRain, chFreqRain,
            chartSDomain, chartZDomain;
    }
}
