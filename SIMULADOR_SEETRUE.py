import customtkinter as ctk
import tkinter as tk
from tkinter import filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Circle
import sys
import os
import math

# --- CONFIGURACIÓN INICIAL ---
ctk.set_appearance_mode("Light")
ctk.set_default_color_theme("blue")

try:
    import rasterio
    from shapely.geometry import LineString
    from pyproj import Transformer, Geod
    HAS_GEOSPATIAL = True
except ImportError:
    HAS_GEOSPATIAL = False
    print("ADVERTENCIA: Librerías geoespaciales no encontradas.")

# --- CONSTANTES FÍSICAS ---
H_PLANCK = 6.626e-34  # J·s
C_LIGHT = 3.0e8       # m/s

class OpticalLinkSimulator:
    def __init__(self, root):
        self.root = root
        self.root.title("Simulador ROBLE -> SEETRUE")
        self.root.geometry("1650x950")
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)
        self.current_theme = "Light"

        # --- VARIABLES DE ESTADO ---
        self.source_type = tk.StringVar(value="laser") 
        self.tx_location = "A" 
        self.show_markers = tk.BooleanVar(value=True) 
        
        # --- RETROCUBO ---
        self.use_retro = tk.BooleanVar(value=False)
        self.var_retro_diam = tk.DoubleVar(value=6.0)
        self.var_retro_eff = tk.DoubleVar(value=90.0)
        
        # --- BASE DE DATOS DE HITOS ---
        self.reference_points = [
            {"name": "Caleu", "lat": -32.998, "lon": -71.045},
            {"name": "La Dormida", "lat": -33.063, "lon": -71.102},
            {"name": "Quebrada Alvarado", "lat": -33.090, "lon": -71.200},
            {"name": "Olmué (Ref)", "lat": -33.00, "lon": -71.30}, 
            {"name": "Villa Alemana", "lat": -33.046, "lon": -71.375},
            {"name": "Peñablanca", "lat": -33.036, "lon": -71.350},
            {"name": "Quilpué", "lat": -33.049, "lon": -71.442},
            {"name": "Lago Peñuelas (Ref)", "lat": -33.15, "lon": -71.48}
        ]

        # Coordenadas por defecto
        self.default_lat_A = -32.976308; self.default_lon_A = -71.014799
        self.default_lat_B = -33.146199; self.default_lon_B = -71.572729
        self.var_name_a = tk.StringVar(value="El Roble")
        self.var_name_b = tk.StringVar(value="Seetrue")

        # Variables Numéricas
        self.var_lat_a = tk.DoubleVar(value=self.default_lat_A)
        self.var_lon_a = tk.DoubleVar(value=self.default_lon_A)
        self.var_lat_b = tk.DoubleVar(value=self.default_lat_B)
        self.var_lon_b = tk.DoubleVar(value=self.default_lon_B)
        self.var_alt_a = tk.DoubleVar(value=2203.0) 
        self.var_alt_b = tk.DoubleVar(value=300.0)  
        
        # Datos del Perfil
        self.dem_path = None
        self.real_profile_dist = None 
        self.real_profile_elev = None 
        self.use_real_profile = False

        # Simulación
        self.var_power_tx = tk.DoubleVar(value=1.5)
        self.var_dist = tk.DoubleVar(value=55.0)
        self.var_diam_tx = tk.DoubleVar(value=2.5)
        self.var_div = tk.DoubleVar(value=5.0)
        self.var_alpha = tk.DoubleVar(value=0.15)
        self.var_lambda = tk.DoubleVar(value=850.0)
        self.var_log_cn2 = tk.DoubleVar(value=-16.0)
        
        # Receptor & Cámara
        self.var_diam_rx = tk.DoubleVar(value=35.0)
        self.var_eff_opt = tk.DoubleVar(value=80.0)
        self.var_qe = tk.DoubleVar(value=70.0)
        self.var_exposure = tk.DoubleVar(value=100.0)
        
        # --- RUIDO ELECTRÓNICO ---
        self.var_read_noise = tk.DoubleVar(value=1.5)
        self.var_dark_current = tk.DoubleVar(value=20.0)

        # Textos informativos
        self.info_diff_limit = tk.StringVar(value="---")
        self.info_direction = tk.StringVar(value="Dirección: A -> B")
        self.info_turb = tk.StringVar(value="---") 
        self.info_rx_location = tk.StringVar(value="3. Medio y Receptor")

        self.setup_ui()
        self.on_source_change() 
        self.update_simulation()

    def setup_ui(self):
        self.main_frame = ctk.CTkFrame(self.root)
        self.main_frame.pack(fill=tk.BOTH, expand=True)

        self.frame_controls = ctk.CTkScrollableFrame(self.main_frame, width=540, label_text="Panel de Control")
        self.frame_controls.pack(side=tk.LEFT, fill=tk.BOTH, padx=10, pady=10)

        self.btn_theme = ctk.CTkButton(self.frame_controls, text="Modo: Día ☀️", 
                                       fg_color="#00BCD4", text_color="black", hover_color="#00ACC1",
                                       command=self.toggle_theme)
        self.btn_theme.pack(fill='x', pady=5)

        self.frame_plot = ctk.CTkFrame(self.main_frame)
        self.frame_plot.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

        # 1. UBICACIÓN
        self.create_section_label("1. Topografía y Ubicación")
        frm_geo = self.create_group_frame("Coordenadas")
        
        f_a = ctk.CTkFrame(frm_geo, fg_color="transparent"); f_a.pack(fill='x')
        ctk.CTkLabel(f_a, text="Lugar A:", width=60).pack(side='left')
        e_a = ctk.CTkEntry(f_a, textvariable=self.var_name_a); e_a.pack(side='left', fill='x', expand=True)
        e_a.bind('<Return>', lambda e: self.update_simulation()); e_a.bind('<FocusOut>', lambda e: self.update_simulation())
        self.create_coord_entry(f_a, "Lat:", self.var_lat_a); self.create_coord_entry(f_a, "Lon:", self.var_lon_a)

        f_b = ctk.CTkFrame(frm_geo, fg_color="transparent"); f_b.pack(fill='x', pady=5)
        ctk.CTkLabel(f_b, text="Lugar B:", width=60).pack(side='left')
        e_b = ctk.CTkEntry(f_b, textvariable=self.var_name_b); e_b.pack(side='left', fill='x', expand=True)
        e_b.bind('<Return>', lambda e: self.update_simulation()); e_b.bind('<FocusOut>', lambda e: self.update_simulation())
        self.create_coord_entry(f_b, "Lat:", self.var_lat_b); self.create_coord_entry(f_b, "Lon:", self.var_lon_b)

        f_btns = ctk.CTkFrame(frm_geo, fg_color="transparent"); f_btns.pack(fill='x', pady=5)
        ctk.CTkButton(f_btns, text="📍 El Roble ↔ Seetrue", command=self.set_preset_roble_seetrue, width=120).pack(side='left', padx=2)
        ctk.CTkButton(f_btns, text="📍 Seetrue ↔ Co. Frutilla", command=self.set_preset_frutilla, width=120).pack(side='left', padx=2)

        frm_dem = ctk.CTkFrame(self.frame_controls, fg_color="transparent"); frm_dem.pack(fill='x', pady=5)
        ctk.CTkButton(frm_dem, text="📂 Cargar DEM (.tif)", command=self.load_dem_file, fg_color="#E0E0E0", text_color="black").pack(fill='x')
        self.lbl_dem_status = ctk.CTkLabel(frm_dem, text="Estado: Terreno sintético", text_color="orange", font=('Arial', 10)); self.lbl_dem_status.pack(anchor='w')
        ctk.CTkCheckBox(frm_dem, text="Mostrar Hitos", variable=self.show_markers, command=self.update_simulation).pack(anchor='w')

        frm_dir = self.create_group_frame("Dirección del Enlace")
        ctk.CTkLabel(frm_dir, textvariable=self.info_direction, font=('Arial', 12, 'bold'), text_color="#1E90FF").pack(anchor='center', pady=5)
        ctk.CTkButton(frm_dir, text="⇄ Invertir Sentido Tx", command=self.toggle_direction, height=24).pack(fill='x')

        # 2. FUENTE
        self.create_section_label("2. Fuente Láser/LED")
        frm_src = ctk.CTkFrame(self.frame_controls, fg_color="transparent"); frm_src.pack(fill='x')
        ctk.CTkRadioButton(frm_src, text="Láser", variable=self.source_type, value="laser", command=self.on_source_change).pack(side='left', padx=10)
        ctk.CTkRadioButton(frm_src, text="LED", variable=self.source_type, value="led", command=self.on_source_change).pack(side='left', padx=10)

        self.create_input_group(self.frame_controls, "Potencia Tx (Watts)", self.var_power_tx, 0.001, 50.0)
        self.create_input_group(self.frame_controls, "Longitud Onda (nm)", self.var_lambda, 380, 1600, step=1.0)
        self.create_input_group(self.frame_controls, "Diámetro Salida (cm)", self.var_diam_tx, 0.1, 50.0)
        ctk.CTkLabel(self.frame_controls, textvariable=self.info_diff_limit, text_color="gray", font=('Arial', 10, 'italic')).pack(anchor='w', padx=5)
        self.create_input_group(self.frame_controls, "Divergencia (mrad)", self.var_div, 0.1, 10.0)

        # 3. RECEPTOR
        ctk.CTkLabel(self.frame_controls, textvariable=self.info_rx_location, font=('Arial', 14, 'bold')).pack(pady=(15,5), anchor='w')
        self.create_input_group(self.frame_controls, "Coef. Extinción (1/km)", self.var_alpha, 0.0, 3.0, step=0.01)
        frm_turb = self.create_group_frame("Turbulencia (log Cn²)")
        self.create_input_group_cn2(frm_turb, self.var_log_cn2, -17.0, -12.0)
        ctk.CTkLabel(frm_turb, textvariable=self.info_turb, font=('Arial', 10, 'bold'), text_color="#d9534f").pack(anchor='w')

        self.create_input_group(self.frame_controls, "Diámetro Cámara (cm)", self.var_diam_rx, 5, 100)
        
        frm_sens = self.create_group_frame("Sensor & Ruido Electrónico")
        self.create_input_group(frm_sens, "Eficiencia Óptica (%)", self.var_eff_opt, 1, 100)
        self.create_input_group(frm_sens, "Eficiencia Cuántica (%)", self.var_qe, 1, 100)
        self.create_input_group(frm_sens, "Exposición (ms)", self.var_exposure, 1, 5000, step=10)
        ctk.CTkLabel(frm_sens, text="--- Componentes de Ruido ---", font=('Arial', 9, 'bold'), text_color="gray").pack(pady=2)
        self.create_input_group(frm_sens, "Ruido Lectura (e- rms)", self.var_read_noise, 0, 10, step=0.1)
        self.create_input_group(frm_sens, "Corriente de Oscuridad (e-/s)", self.var_dark_current, 0, 1000, step=1.0) 
        
        # 4. RETROCUBO
        self.create_section_label("4. Retrocubo")
        frm_retro = ctk.CTkFrame(self.frame_controls, border_width=2, border_color="#FFD700") 
        frm_retro.pack(fill='x', pady=5, padx=2)
        sw_retro = ctk.CTkSwitch(frm_retro, text="ACTIVAR Retrocubo", variable=self.use_retro, 
                                 command=self.update_simulation, font=('Arial', 12, 'bold'))
        sw_retro.pack(anchor='w', padx=10, pady=10)
        self.create_input_group(frm_retro, "Diámetro Retro (cm)", self.var_retro_diam, 1.0, 50.0)
        self.create_input_group(frm_retro, "Reflectancia (%)", self.var_retro_eff, 1.0, 100.0)
        ctk.CTkLabel(frm_retro, text="Nota: Rx pasa al origen.", font=('Arial', 10, 'italic')).pack(pady=5)

        # RESULTADOS
        self.result_frame = self.create_group_frame("Resultados")
        self.lbl_res_power = ctk.CTkLabel(self.result_frame, text="---", font=('Courier', 16, 'bold'), text_color="#1E90FF")
        self.lbl_res_power.pack(anchor="w", pady=5)
        self.lbl_scint = ctk.CTkLabel(self.result_frame, text="---")
        self.lbl_scint.pack(anchor="w")
        self.lbl_irrad = ctk.CTkLabel(self.result_frame, text="---")
        self.lbl_irrad.pack(anchor="w")
        self.lbl_spot = ctk.CTkLabel(self.result_frame, text="---")
        self.lbl_spot.pack(anchor="w")
        
        self.lbl_noise_breakdown = ctk.CTkLabel(self.result_frame, text="", text_color="gray", font=('Arial', 10))
        self.lbl_noise_breakdown.pack(anchor="w")
        self.lbl_snr = ctk.CTkLabel(self.result_frame, text="SNR: ---", font=('Arial', 14, 'bold'))
        self.lbl_snr.pack(anchor="w", pady=5)
        self.lbl_retro_info = ctk.CTkLabel(self.result_frame, text="", text_color="#FF8C00", font=('Arial', 11))
        self.lbl_retro_info.pack(anchor="w")

        self.fig, (self.ax_terrain, self.ax_spot) = plt.subplots(2, 1, figsize=(8, 9), gridspec_kw={'height_ratios': [2, 1]})
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame_plot)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # --- HELPERS UI ---
    def create_section_label(self, text):
        ctk.CTkLabel(self.frame_controls, text=text, font=('Arial', 14, 'bold')).pack(pady=(15,5), anchor='w')
    def create_group_frame(self, title):
        frm = ctk.CTkFrame(self.frame_controls)
        frm.pack(fill='x', pady=5, padx=2)
        ctk.CTkLabel(frm, text=title, font=('Arial', 11, 'bold'), text_color="gray").pack(anchor='w', padx=5, pady=(5,0))
        return frm
    def create_coord_entry(self, parent, label, variable):
        ctk.CTkLabel(parent, text=label).pack(side='left', padx=(5,2))
        e = ctk.CTkEntry(parent, textvariable=variable, width=80); e.pack(side='left')
        e.bind('<Return>', lambda e: self.process_geography())
    def create_input_group(self, parent, label_text, variable, min_val, max_val, step=0.1):
        frame = ctk.CTkFrame(parent, fg_color="transparent"); frame.pack(fill="x", pady=2)
        var_str = tk.StringVar(value=f"{variable.get():g}")
        top = ctk.CTkFrame(frame, fg_color="transparent"); top.pack(fill="x")
        ctk.CTkLabel(top, text=label_text, font=('Arial', 11)).pack(side="left")
        def on_entry_commit(event=None):
            try: val = float(var_str.get()); variable.set(val); self.update_simulation()
            except: pass
        entry = ctk.CTkEntry(top, textvariable=var_str, width=60); entry.pack(side="right")
        entry.bind('<Return>', on_entry_commit); entry.bind('<FocusOut>', on_entry_commit)
        steps = int((max_val - min_val) / step)
        scale = ctk.CTkSlider(frame, from_=min_val, to=max_val, number_of_steps=steps,
                              command=lambda v: [variable.set(float(v)), self.update_simulation()])
        scale.set(variable.get()); scale.pack(fill="x", pady=(2,0))
        variable.trace_add("write", lambda *a: [scale.set(variable.get()), var_str.set(f"{variable.get():g}") if self.root.focus_get() != entry._entry else None])
        return [entry, scale]
    def create_input_group_cn2(self, parent, variable, min_val, max_val):
        self.create_input_group(parent, "Valor Log (10^x):", variable, min_val, max_val, step=0.1)

    # --- SIMULACIÓN FÍSICA RIGUROSA ---
    def update_simulation(self):
        try:
            wl_nm = self.var_lambda.get(); lam_m = wl_nm * 1e-9; k = 2 * np.pi / lam_m
            L_km = self.var_dist.get(); L_m = L_km * 1000.0
            Cn2_ground = 10**(self.var_log_cn2.get()); alpha = self.var_alpha.get()
            name_A = self.var_name_a.get(); name_B = self.var_name_b.get()
            
            if self.tx_location == "A": txt_dir = f"Tx: {name_A} ➔ Rx: {name_B}"
            else: txt_dir = f"Tx: {name_B} ➔ Rx: {name_A}"
            self.info_direction.set(f"{txt_dir}\nDistancia Óptica: {L_km:.3f} km")
            
            # Turbulencia Base (Hufnagel simple)
            alt_A = self.var_alt_a.get(); alt_B = self.var_alt_b.get()
            h_avg = max((alt_A + alt_B)/2.0, 10.0)
            Cn2_eff = Cn2_ground * np.exp(-h_avg / 3000.0)
            sigma_R2 = 1.23 * Cn2_eff * (k**(7/6)) * (L_m**(11/6)) # Variancia de Rytov 
            
            turb_str = "Débil" if sigma_R2 < 1 else "Fuerte"
            self.info_turb.set(f"Rytov σ_R²: {sigma_R2:.2f} ({turb_str})")

            # --- IDA (Transmisor -> Destino) ---
            P_tx = self.var_power_tx.get(); D_tx_m = self.var_diam_tx.get() / 100.0
            if self.source_type.get() == "laser":
                theta_diff = (2 * lam_m) / (np.pi * D_tx_m)
                theta_geo = (self.var_div.get() / 1000.0) / 2.0
                theta_tot_fwd = np.sqrt(theta_diff**2 + theta_geo**2)
                self.info_diff_limit.set(f"Difracción: {(theta_diff*1000):.3f} mrad")
            else:
                theta_tot_fwd = np.radians(self.var_div.get()) / 2.0
                self.info_diff_limit.set("N/A (LED)")

            w0_fwd = D_tx_m / 2.0
            W_vac_fwd = w0_fwd + L_m * np.tan(theta_tot_fwd)
            W_LT_fwd = W_vac_fwd * np.sqrt(1 + 0.5*sigma_R2)
            Area_beam_fwd = np.pi * (W_LT_fwd**2)
            trans_atm = np.exp(-alpha * L_km)
            
            if self.source_type.get() == "led": Irr_peak_fwd = (P_tx * trans_atm) / Area_beam_fwd
            else: Irr_peak_fwd = (2 * P_tx * trans_atm) / Area_beam_fwd

            # --- PARÁMETROS FINALES Y RETROCUBO ---
            D_rx_m = self.var_diam_rx.get() / 100.0
            
            if self.use_retro.get():
                self.info_rx_location.set(f"3. Cámara en {name_A if self.tx_location == 'A' else name_B} (Retorno)")
                D_retro_m = self.var_retro_diam.get() / 100.0
                Area_retro = np.pi * ((D_retro_m/2)**2)
                eff_retro = self.var_retro_eff.get() / 100.0
                
                # Potencia capturada (Ida)
                P_arriving_retro = min(Irr_peak_fwd * Area_retro, P_tx * trans_atm)
                P_reflected = P_arriving_retro * eff_retro
                self.lbl_retro_info.configure(text=f"Al Retro: {self.format_power(P_arriving_retro)[1]:.2f} {self.format_power(P_arriving_retro)[0]} | Reflejado: {self.format_power(P_reflected)[1]:.2f} {self.format_power(P_reflected)[0]}")

                # Vuelta (Difracción de Retro)
                theta_ret = (4 * lam_m) / (np.pi * D_retro_m)
                w0_ret = D_retro_m / 2.0
                W_vac_ret = w0_ret + L_m * np.tan(theta_ret)
                W_LT_ret = W_vac_ret * np.sqrt(1 + 0.5*sigma_R2) # Turbulencia afecta radio también
                Area_beam_ret = np.pi * (W_LT_ret**2)
                Irr_peak_ret = (2 * P_reflected * trans_atm) / Area_beam_ret 
                
                Final_Irradiance = Irr_peak_ret
                Final_W_LT = W_LT_ret
                plot_W_end = W_LT_fwd; plot_W_ret_start = D_retro_m/2.0; plot_W_ret_end = W_LT_ret
                
                # === Física de Retrocubo
                factor_A_fwd = 1 / (1 + (np.pi * D_retro_m**2)/(4*lam_m*L_m))
                sigma_I2_fwd = min(sigma_R2 * factor_A_fwd, 1.2)
                
                # 2. Centelleo Vuelta (Visto por el Telescopio)
                # Apertura = D_rx 
                factor_A_ret = 1 / (1 + (np.pi * D_rx_m**2)/(4*lam_m*L_m))
                sigma_I2_ret = min(sigma_R2 * factor_A_ret, 1.2)
                
                # 3. Centelleo Total (Varianzas Independientes)
                sigma_I2_total = (1 + sigma_I2_fwd)*(1 + sigma_I2_ret) - 1
                
            else:
                self.info_rx_location.set(f"3. Receptor en {name_B if self.tx_location == 'A' else name_A} (Destino)")
                self.lbl_retro_info.configure(text="")
                Final_Irradiance = Irr_peak_fwd
                Final_W_LT = W_LT_fwd
                plot_W_end = W_LT_fwd; plot_W_ret_start = None; plot_W_ret_end = None
                
                # Centelleo Estándar
                factor_A = 1 / (1 + (np.pi * D_rx_m**2)/(4*lam_m*L_m))
                sigma_I2_total = min(sigma_R2 * factor_A, 1.2)

            # === DETECCIÓN y SNR ===
            Area_rx = np.pi * ((D_rx_m/2)**2)
            eff_opt = self.var_eff_opt.get() / 100.0
            P_rx_mean = Final_Irradiance * Area_rx * eff_opt
            
            E_ph = (H_PLANCK * C_LIGHT) / lam_m
            flux_ph = P_rx_mean / E_ph
            t_exp_ms = self.var_exposure.get(); t_exp_s = t_exp_ms / 1000.0
            signal = flux_ph * t_exp_s * (self.var_qe.get()/100.0)
            
            #== Ruido
            sigma_read = self.var_read_noise.get()
            var_read = sigma_read**2
            dark_rate = self.var_dark_current.get()
            var_dark = dark_rate * t_exp_s
            var_shot = signal
            
            #== Centelleo
            var_scint = (signal * np.sqrt(sigma_I2_total))**2
            
            total_var = var_shot + var_read + var_dark + var_scint
            noise = np.sqrt(total_var)
            snr = signal / noise if noise > 0 else 0

            if snr > 0:
                cv_percent = (1.0 / snr) * 100.0
                ber_val = 0.5 * math.erfc(snr / np.sqrt(2))
            else:
                cv_percent = 100.0 
                ber_val = 0.5      

            # UI Resultados (Actualización de etiquetas)
            unit, val = self.format_power(P_rx_mean)
            txt_modo = "(RETORNO)" if self.use_retro.get() else "(IDA)"
            self.lbl_res_power.configure(text=f"Potencia Rx {txt_modo}: {val:.2f} {unit}")
            self.lbl_irrad.configure(text=f"Irradiancia: {Final_Irradiance:.2e} W/m²")
            self.lbl_spot.configure(text=f"Spot Final (LT): {Final_W_LT*2:.1f} m")
            self.lbl_scint.configure(text=f"Centelleo Total σ_I²: {sigma_I2_total:.3f}")
            
            # Desglose de Ruido
            self.lbl_noise_breakdown.configure(text=f"Ruido(e-): Read={sigma_read:.1f} | Dark={np.sqrt(var_dark):.1f} | Shot={np.sqrt(var_shot):.1f} | Scint={np.sqrt(var_scint):.1f}")
            
            # Colores del semáforo basados en BER/SNR
            col = "#00FF00" if snr > 10 else "#FFA500" if snr > 3 else "#FF4444"
            
            # --- LÍNEA FINAL MODIFICADA ---
            # Muestra SNR, Porcentaje de Ruido (CV) y Tasa de Error (BER) en notación científica
            self.lbl_snr.configure(text=f"SNR: {snr:.2f}  |  CV: {cv_percent:.1f}%  |  BER: {ber_val:.1e}", text_color=col)
            
            # Plotting
            lat_A, lon_A = self.var_lat_a.get(), self.var_lon_a.get()
            lat_B, lon_B = self.var_lat_b.get(), self.var_lon_b.get()
            markers = self.get_markers_on_path(lat_A, lon_A, lat_B, lon_B, L_km)
            
            self.plot_terrain_v16(L_km, plot_W_end, wl_nm, alt_A, alt_B, D_tx_m, markers, plot_W_ret_start, plot_W_ret_end, name_A, name_B)
            self.plot_spot(Final_W_LT, D_rx_m, wl_nm)
            self.canvas.draw()
        except Exception as e: print(e)

    # --- PLOTTING HELPERS ---
    def plot_terrain_v16(self, L_km, W_radius_end, wl_nm, alt_A, alt_B, D_tx_m, markers, ret_start=None, ret_end=None, name_A="A", name_B="B"):
        ax = self.ax_terrain; ax.clear()
        c_beam = self.wavelength_to_hex(wl_nm)
        if self.use_real_profile and self.real_profile_dist is not None:
            x_terr = self.real_profile_dist; y_terr = self.real_profile_elev
            min_elev = np.nanmin(y_terr); max_elev = np.nanmax(y_terr)
            ax.fill_between(x_terr, 0, y_terr, color='#4E342E', alpha=0.9)
        else:
            x_terr = np.linspace(0, L_km, 200)
            base_line = alt_A + (alt_B - alt_A)*(x_terr/L_km)
            y_terr = base_line - 800 * np.sin(np.pi * x_terr/L_km)**2
            y_terr = np.clip(y_terr, 0, None)
            min_elev = min(alt_A, alt_B); max_elev = max(alt_A, alt_B)
            ax.fill_between(x_terr, 0, y_terr, color='#5D4037', alpha=0.7)

        x_beam = np.linspace(0, L_km, 100); los = alt_A + (alt_B - alt_A)*(x_beam/L_km); R_start = D_tx_m / 2.0
        if self.tx_location == "A": prof = R_start + (W_radius_end - R_start) * (x_beam / L_km); lbl_ida = f"Ida ({name_A} -> {name_B})"
        else: prof = W_radius_end + (R_start - W_radius_end) * (x_beam / L_km); lbl_ida = f"Ida ({name_B} -> {name_A})"
        ax.fill_between(x_beam, los-prof, los+prof, color=c_beam, alpha=0.4, label=lbl_ida)

        if self.use_retro.get() and ret_start is not None:
            c_ret = "#00FFFF" 
            if self.tx_location == "A": prof_ret = ret_end + (ret_start - ret_end) * (x_beam / L_km)
            else: prof_ret = ret_start + (ret_end - ret_start) * (x_beam / L_km)
            ax.plot(x_beam, los+prof_ret, color=c_ret, linestyle='--', linewidth=1, alpha=0.8)
            ax.plot(x_beam, los-prof_ret, color=c_ret, linestyle='--', linewidth=1, alpha=0.8)
            ax.fill_between(x_beam, los-prof_ret, los+prof_ret, color=c_ret, alpha=0.2, label='Haz Retorno')

        ax.annotate(f"{name_A}\n{alt_A:.0f}m", xy=(0, alt_A), xytext=(0, 20), textcoords='offset points', ha='center', fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.7))
        ax.annotate(f"{name_B}\n{alt_B:.0f}m", xy=(L_km, alt_B), xytext=(0, 20), textcoords='offset points', ha='center', fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.7))
        if markers:
            for d, n in markers:
                ax.vlines(d, 0, max_elev, colors='gray', linestyles=':', alpha=0.5); ax.text(d, max_elev+100, n, rotation=90, fontsize=8)
        ax.set_ylim(max(0, min_elev - 200), max_elev + 1200)
        ax.set_title(f"Perfil del Enlace {'(Con Rebote)' if self.use_retro.get() else ''}")
        ax.set_ylabel("Altitud (m)"); ax.set_xlabel("Distancia (km)"); ax.legend(loc='upper right', fontsize='8'); ax.grid(True, alpha=0.3)

    def plot_spot(self, W_radius, D_rx, wl_nm):
        ax = self.ax_spot; ax.clear()
        c_beam = "#00FFFF" if self.use_retro.get() else self.wavelength_to_hex(wl_nm)
        lbl_beam = 'Haz Retorno' if self.use_retro.get() else 'Haz Incidente'
        beam = Circle((0,0), W_radius, color=c_beam, alpha=0.5, label=lbl_beam)
        tel = Circle((0,0), D_rx/2, edgecolor='black', facecolor='none', lw=2, ls='--', label='Cámara Rx')
        ax.add_patch(beam); ax.add_patch(tel)
        limit = max(D_rx*2, min(W_radius*1.2, W_radius))
        ax.set_xlim(-limit, limit); ax.set_ylim(-limit, limit)
        ax.set_aspect('equal'); ax.set_title("Plano del Detector")
        ax.set_xlabel("Posición (m)"); ax.set_ylabel("Posición (m)")
        
        # Calcular Areas
        area_spot = np.pi * W_radius**2
        area_rx = np.pi * (D_rx/200.0)**2
        ratio = area_rx / area_spot if area_spot > 0 else 0
        txt_areas = f"Area Spot: {area_spot:.1f} m² | Area Rx: {area_rx*1e4:.1f} cm² | Captura: {ratio:.2e}"
        ax.text(0.5, -0.15, txt_areas, transform=ax.transAxes, ha='center', fontsize=9, bbox=dict(facecolor='white', alpha=0.8))
        
        ax.grid(True, linestyle='--', alpha=0.3); ax.legend(fontsize='small', loc='upper right')

    def format_power(self, p):
        if p > 1e-3: return "mW", p*1e3
        if p > 1e-6: return "µW", p*1e6
        if p > 1e-9: return "nW", p*1e9
        if p > 1e-12: return "pW", p*1e12
        if p > 1e-15: return "fW", p*1e15
        if p > 0: return "aW", p*1e18
        return "pW", 0.0

    def wavelength_to_hex(self, wl_nm):
        wl = float(wl_nm)
        if wl >= 380 and wl <= 780:
            if wl < 440: r,g,b = -(wl-440)/(440-380), 0, 1
            elif wl < 490: r,g,b = 0, (wl-440)/(490-440), 1
            elif wl < 510: r,g,b = 0, 1, -(wl-510)/(510-490)
            elif wl < 580: r,g,b = (wl-510)/(580-510), 1, 0
            elif wl < 645: r,g,b = 1, -(wl-645)/(645-580), 0
            else: r,g,b = 1, 0, 0
            return f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}'
        elif wl > 780: return "#8B0000"
        else: return "#4B0082"
    def on_source_change(self): self.update_simulation()
    def toggle_theme(self):
        self.current_theme = "Dark" if self.current_theme == "Light" else "Light"
        ctk.set_appearance_mode(self.current_theme)
        self.btn_theme.configure(text=f"Modo: {'Noche 🌙' if self.current_theme == 'Dark' else 'Día ☀️'}")
    def toggle_direction(self):
        self.tx_location = "B" if self.tx_location == "A" else "A"
        self.update_simulation()
    def load_dem_file(self):
        if not HAS_GEOSPATIAL: return
        fn = filedialog.askopenfilename(filetypes=[("GeoTIFF", "*.tif")])
        if fn:
            self.dem_path = fn
            self.lbl_dem_status.configure(text=f"Archivo: {os.path.basename(fn)}", text_color="green")
            self.process_geography()
    def process_geography(self):
        lat1, lon1 = self.var_lat_a.get(), self.var_lon_a.get()
        lat2, lon2 = self.var_lat_b.get(), self.var_lon_b.get()
        dist_km = 55.0
        if HAS_GEOSPATIAL:
            try:
                geod = Geod(ellps="WGS84"); _, _, total_dist_m = geod.inv(lon1, lat1, lon2, lat2); dist_km = total_dist_m / 1000.0
            except: pass
        self.var_dist.set(dist_km)
        if HAS_GEOSPATIAL and self.dem_path:
            try:
                with rasterio.open(self.dem_path) as src:
                    transformer = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True)
                    x1, y1 = transformer.transform(lon1, lat1); x2, y2 = transformer.transform(lon2, lat2)
                    line = LineString([(x1, y1), (x2, y2)]); dists = np.linspace(0, line.length, 500)
                    points = [line.interpolate(d) for d in dists]; coords = [(p.x, p.y) for p in points]
                    samples = list(src.sample(coords))
                    elevs = np.array([s[0] if s[0] != src.nodata else np.nan for s in samples])
                    nans, x = np.isnan(elevs), lambda z: z.nonzero()[0]
                    if np.any(nans): elevs[nans] = np.interp(x(nans), x(~nans), elevs[~nans])
                    self.real_profile_elev = elevs; self.real_profile_dist = np.linspace(0, dist_km, 500)
                    self.use_real_profile = True; self.var_alt_a.set(elevs[0]); self.var_alt_b.set(elevs[-1])
            except: self.use_real_profile = False
        else: self.use_real_profile = False
        self.update_simulation()
    def get_markers_on_path(self, lat_A, lon_A, lat_B, lon_B, total_km):
        if not HAS_GEOSPATIAL or not self.show_markers.get(): return []
        found = []; geod = Geod(ellps="WGS84")
        for m in self.reference_points:
            _, _, d_A = geod.inv(lon_A, lat_A, m["lon"], m["lat"]); d_km = d_A / 1000.0
            _, _, d_B = geod.inv(m["lon"], m["lat"], lon_B, lat_B); detour = (d_A + d_B)/1000.0
            if abs(detour - total_km) < 15.0 and 0.5 < d_km < (total_km - 0.5): found.append((d_km, m["name"]))
        return found
    def set_preset_roble_seetrue(self):
        self.var_lat_a.set(self.default_lat_A); self.var_lon_a.set(self.default_lon_A)
        self.var_lat_b.set(self.default_lat_B); self.var_lon_b.set(self.default_lon_B)
        self.var_name_a.set("El Roble"); self.var_name_b.set("Seetrue")
        self.process_geography()
    def set_preset_frutilla(self):
        self.var_lat_a.set(-33.161667); self.var_lon_a.set(-71.432500)
        self.var_lat_b.set(self.default_lat_B); self.var_lon_b.set(self.default_lon_B)
        self.var_name_a.set("Cerro Frutilla"); self.var_name_b.set("Seetrue")
        self.process_geography()
    def on_close(self): self.root.quit(); self.root.destroy(); sys.exit()

if __name__ == "__main__":
    root = ctk.CTk()
    app = OpticalLinkSimulator(root)
    root.mainloop()