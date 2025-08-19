#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT
import pandas as pd
import numpy as np
import mplhep as mlp
import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import math

from ROOT import (
    RooRealVar, RooDataHist, RooArgList,
    RooCBShape, RooPolynomial, RooAddPdf, RooFit,
    TPaveText, TFile, TCanvas,  RooGaussian
)


from ROOT import TPaveText, gStyle


from ROOT import (
    RooRealVar, RooDataHist, RooArgList, RooCBShape, RooVoigtian, RooAddPdf,
    RooFit, TCanvas, TPaveText, TLegend, kBlack, kRed, kBlue, kDashed
)


# In[2]:


ROOT.ROOT.EnableImplicitMT()

file_path = "root://eosuser.cern.ch//eos/user/j/jobarrei/realdata/dados_cbpf_hhh/*.root"
tree_name = "Turbo_BuToKpKpKm/DecayTree"

df_kkk = ROOT.RDataFrame(tree_name, file_path)
df_kkpi = ROOT.RDataFrame(tree_name, file_path)
df_pipipi = ROOT.RDataFrame(tree_name, file_path)


# In[3]:


df_kkk = (
    df_kkk
    
    # Corte de vertex
    #.Filter("Bp_BPVIPCHI2 < 10")
    
    
    # =========================
    # Cortes completos para K+ (Kp1)
    # =========================
    .Filter("Kp1_PT > 500")
    .Filter("Kp1_ETA > 2 && Kp1_ETA < 5")
    .Filter("Kp1_P > 1000 && Kp1_P < 100000")
    .Filter("Kp1_TRCHI2DOF < 3")
    .Filter("Kp1_TRGHOSTPROB < 0.5")
    .Filter("Kp1_ISMUON == 0")
    
    
    # =========================
    # Cortes completos para K+ (Kp2)
    # =========================
    .Filter("Kp2_PT > 500")
    .Filter("Kp2_ETA > 2 && Kp2_ETA < 5")
    .Filter("Kp2_P > 1000 && Kp2_P < 100000")
    .Filter("Kp2_TRCHI2DOF < 3")
    .Filter("Kp2_TRGHOSTPROB < 0.5")
    .Filter("Kp2_ISMUON == 0")
    
    # =========================
    # Cortes completos para K- (Km)
    # =========================
    
    .Filter("Km_PT > 500")
    .Filter("Km_ETA > 2 && Km_ETA < 5")
    .Filter("Km_P > 1000 && Km_P < 100000")
    .Filter("Km_TRCHI2DOF < 3")
    .Filter("Km_TRGHOSTPROB < 0.5")
    .Filter("Km_ISMUON == 0")
    
    
    
    #Cortes no PID
    .Filter("Kp1_PID_K > 5")
    .Filter("Kp2_PID_K > 5")
    .Filter("Km_PID_K > 5")
    
    
    
    # =========================
    # Cortes gerais do B+
    # =========================
    
    .Filter("Bp_BPVFDCHI2 > 500")
    .Filter("Bp_VCHI2DOF < 5")
    .Filter("Bp_BPVIPCHI2 < 10")
    
    # Cortes de massa e momentum
    .Filter("Bp_PT > 1000")
    .Filter("Bp_M > 5050 && Bp_M < 6300")
    .Filter("Bp_MCORR > 4000 && Bp_MCORR < 7000")
    
    # Corte de direção e geometria
    .Filter("Bp_BPVDIRA > 0.99998")
    .Filter("Bp_MAXDOCA < 0.2")
    
    # Cortes combinados de momentum das tracks
    .Filter("Bp_SUMPT > 4500")
    .Filter("Bp_SUMP > 20000")
    .Filter("Bp_MAXPT > 1500")
    
    # Definir coluna de massa padronizada (opcional)
    .Define("Bp_m", "Bp_M")
)


# In[4]:


c1 = ROOT.TCanvas("c1", "Canvas", 800, 600)


h_Bp_mass_kkk = df_kkk.Histo1D(
    ("Bp_mass", "Distribuição de Massa do B^{+}; M_{K^{+}K^{-}K^{+}} [MeV]; Eventos", 100, 5150, 5400),
    "Bp_m"
)

h_Bp_mass_kkk.Draw("HIST")
c1.Update()
c1.Draw()


# In[5]:


from ROOT import (
    TFile, RooRealVar, RooDataHist, RooCBShape, RooExponential, RooAddPdf,
    RooArgList, RooArgSet, RooFit, TCanvas, TLatex, TLegend, kRed, kBlue, kGreen, kBlack
)

# --- Abrir histograma ---
f = TFile.Open("B_KKK_hist.root")
hist = f.Get("h_Bp_mass_kkk")

# --- Variável de massa em MeV ---
x = RooRealVar("x", "Massa [MeV]", 5150, 5400)
x.setBins(hist.GetNbinsX())

# --- Dataset RooFit ---
data = RooDataHist("data", "Dados do histograma", RooArgList(x), RooFit.Import(hist))

# --- Crystal Balls ---
m_cb1 = RooRealVar("m_cb1", "Massa CB1", 5279, 5250, 5300)
sigma_cb1 = RooRealVar("sigma_cb1", "Largura CB1", 10, 5, 20)
alpha1 = RooRealVar("alpha1", "alpha1", 1.5, 1.0, 3.0)
n1 = RooRealVar("n1", "n1", 2.0, 1.0, 10.0)
cb_pdf1 = RooCBShape("cb_pdf1", "Crystal Ball 1", x, m_cb1, sigma_cb1, alpha1, n1)

m_cb2 = RooRealVar("m_cb2", "Massa CB2", 5270, 5250, 5300)
sigma_cb2 = RooRealVar("sigma_cb2", "Largura CB2", 15, 5, 20)
alpha2 = RooRealVar("alpha2", "alpha2", 1.5, 1.0, 3.0)
n2 = RooRealVar("n2", "n2", 2.0, 1.0, 10.0)
cb_pdf2 = RooCBShape("cb_pdf2", "Crystal Ball 2", x, m_cb2, sigma_cb2, alpha2, n2)

# --- Fundo exponencial ajustado (mais decaimento) ---
c_bkg = RooRealVar("c_bkg", "Inclinação do fundo", -0.0008, -0.005, -0.0001)
bkg_pdf = RooExponential("bkg_pdf", "Fundo exponencial", x, c_bkg)

# --- Normalizações ---
n_cb1 = RooRealVar("n_cb1", "Eventos CB1", 1000, 100, 10000)
n_cb2 = RooRealVar("n_cb2", "Eventos CB2", 1000, 100, 10000)
n_bkg = RooRealVar("n_bkg", "Eventos de fundo", 500, 100, 5000)

# --- Modelo total ---
model = RooAddPdf("model", "Modelo total",
                  RooArgList(cb_pdf1, cb_pdf2, bkg_pdf),
                  RooArgList(n_cb1, n_cb2, n_bkg))

# --- Fit ---
fit_result = model.fitTo(data, RooFit.PrintLevel(-1), RooFit.Save())

# --- Canvas e frame ---
c1 = TCanvas("c1", "Fit B^{+} -> K^{+}K^{-}K^{+}", 900, 700)
frame = x.frame(RooFit.Title("Distribuição de Massa do B^{+}"))
data.plotOn(frame, RooFit.Name("Dados"), RooFit.MarkerStyle(20), RooFit.MarkerColor(kBlack))
model.plotOn(frame, RooFit.Name("Modelo Total"), RooFit.LineColor(kBlack), RooFit.LineWidth(2))
model.plotOn(frame, RooFit.Components(RooArgSet(cb_pdf1)), RooFit.LineStyle(2), RooFit.LineColor(kRed), RooFit.Name("CB1"))
model.plotOn(frame, RooFit.Components(RooArgSet(cb_pdf2)), RooFit.LineStyle(2), RooFit.LineColor(kBlue), RooFit.Name("CB2"))
model.plotOn(frame, RooFit.Components(RooArgSet(bkg_pdf)), RooFit.LineStyle(2), RooFit.LineColor(kGreen+2), RooFit.Name("Fundo"))

# --- Legenda ---
legend = TLegend(0.65, 0.6, 0.87, 0.87)
legend.SetBorderSize(0)
legend.SetTextSize(0.03)
legend.AddEntry(frame.findObject("Dados"), "Dados", "pe")
legend.AddEntry(frame.findObject("Modelo Total"), "Modelo Total", "l")
legend.AddEntry(frame.findObject("CB1"), "Crystal Ball 1", "l")
legend.AddEntry(frame.findObject("CB2"), "Crystal Ball 2", "l")
legend.AddEntry(frame.findObject("Fundo"), "Fundo Exponencial", "l")

# --- Chi2/NDF ---
chi2 = frame.chiSquare("Modelo Total", "Dados", fit_result.floatParsFinal().getSize())
ndf = hist.GetNbinsX() - fit_result.floatParsFinal().getSize()
chi2_ = chi2/ndf
frame.GetYaxis().SetLabelSize(0.025)

# --- Caixa de parâmetros ---
tex = TLatex()
tex.SetNDC()
tex.SetTextSize(0.03)
tex.SetTextFont(42)

x_tex = 0.6
y_tex = 0.8
dy = 0.04

tex.DrawLatex(x_tex, y_tex, f"CB1: m = {m_cb1.getVal():.1f} #pm {m_cb1.getError():.1f} MeV, #sigma = {sigma_cb1.getVal():.1f} MeV")
tex.DrawLatex(x_tex, y_tex - dy, f"CB2: m = {m_cb2.getVal():.1f} #pm {m_cb2.getError():.1f} MeV, #sigma = {sigma_cb2.getVal():.1f} MeV")
tex.DrawLatex(x_tex, y_tex - 2*dy, f"Eventos: CB1 = {n_cb1.getVal():.0f}, CB2 = {n_cb2.getVal():.0f}, Fundo = {n_bkg.getVal():.0f}")
tex.DrawLatex(x_tex, y_tex - 3*dy, f"Total = {n_cb1.getVal()+n_cb2.getVal()+n_bkg.getVal():.0f}")
tex.DrawLatex(x_tex, y_tex - 4*dy, f"#chi^2/NDF = {chi2:.2f}")

# --- Desenhar ---
frame.Draw()
legend.Draw()
c1.Update()
c1.Draw()


# In[ ]:





# In[6]:


from ROOT import gStyle

# --- Canvas ---
c1 = TCanvas("c1", "Fit B^{+} -> K^{+}K^{-}K^{+}", 900, 700)
c1.cd()

# --- Ajuste do eixo y para múltiplos de 10^3 ---
frame.SetMinimum(0)
frame.GetYaxis().SetLabelSize(0.03)
frame.GetYaxis().SetTitleOffset(1.2)

# --- Desenhar frame e legendas primeiro ---
frame.Draw()
legend.Draw()

# --- Caixa de parâmetros ---
tex = TLatex()
tex.SetNDC()
tex.SetTextSize(0.026)
tex.SetTextFont(40)

frame.GetYaxis().SetLabelSize(0.025)


x_tex = 0.6
y_tex = 0.8
dy = 0.04

tex.DrawLatex(x_tex, y_tex, f"CB1: m = {m_cb1.getVal():.1f} #pm {m_cb1.getError():.1f} MeV, #sigma = {sigma_cb1.getVal():.1f} MeV")
tex.DrawLatex(x_tex, y_tex - dy, f"CB2: m = {m_cb2.getVal():.1f} #pm {m_cb2.getError():.1f} MeV, #sigma = {sigma_cb2.getVal():.1f} MeV")
#tex.DrawLatex(x_tex, y_tex - 2*dy, f"Eventos: CB1 = {n_cb1.getVal():.0f}, CB2 = {n_cb2.getVal():.0f}, Fundo = {n_bkg.getVal():.0f}")
#tex.DrawLatex(x_tex, y_tex - 3*dy, f"Total = {n_cb1.getVal()+n_cb2.getVal()+n_bkg.getVal():.0f}")
tex.DrawLatex(x_tex, y_tex - 4*dy, f"#chi{2}/ndf = {chi2:.2f}")

c1.Update()
c1.Draw()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[7]:


df_kkpi = (
    df_kkpi

    # =========================
    # Cortes completos para pi+ (Kp1)
    # =========================
    .Filter("Kp1_PT > 500")
    .Filter("Kp1_ETA > 2 && Kp1_ETA < 5")
    .Filter("Kp1_P > 1000 && Kp1_P < 100000")
    .Filter("Kp1_TRCHI2DOF < 3")
    .Filter("Kp1_TRGHOSTPROB < 0.5")
    .Filter("Kp1_ISMUON == 0")
    .Filter("Kp1_PID_K < -5")  # garante que não é kaon

    # =========================
    # Cortes completos para K+ (Kp2)
    # =========================
    .Filter("Kp2_PT > 500")
    .Filter("Kp2_ETA > 2 && Kp2_ETA < 5")
    .Filter("Kp2_P > 1000 && Kp2_P < 100000")
    .Filter("Kp2_TRCHI2DOF < 3")
    .Filter("Kp2_TRGHOSTPROB < 0.5")
    .Filter("Kp2_ISMUON == 0")
    .Filter("Kp2_PID_K > 5")     # garante que é kaon

    # =========================
    # Cortes completos para K- (Km)
    # =========================
    .Filter("Km_PT > 500")
    .Filter("Km_ETA > 2 && Km_ETA < 5")
    .Filter("Km_P > 1000 && Km_P < 100000")
    .Filter("Km_TRCHI2DOF < 3")
    .Filter("Km_TRGHOSTPROB < 0.5")
    .Filter("Km_ISMUON == 0")
    .Filter("Km_PID_K > 5")  # garante que é kaon

    # =========================
    # Cortes gerais do B+
    # =========================
    .Filter("Bp_BPVFDCHI2 > 500")
    .Filter("Bp_VCHI2DOF < 5")
    .Filter("Bp_BPVIPCHI2 < 10")
    .Filter("Bp_PT > 1000")
    .Filter("Bp_M > 5050 && Bp_M < 6300")
    .Filter("Bp_MCORR > 4000 && Bp_MCORR < 7000")
    .Filter("Bp_BPVDIRA > 0.99998")
    .Filter("Bp_MAXDOCA < 0.2")
    .Filter("Bp_SUMPT > 4500")
    .Filter("Bp_SUMP > 20000")
    .Filter("Bp_MAXPT > 1500")

    # =========================
    # Definir coluna de massa padronizada
    # =========================
    .Define("Bp_m", "Bp_M")
)


# In[8]:


c2 = ROOT.TCanvas("c2", "Canvas", 800, 600)


h_Bp_mass_kkpi = df_kkpi.Histo1D(
    ("Bp_mass", "Distribuição de Massa do B^{+}; M_{#pi^{+}K^{-}K^{+}} [MeV]; Eventos", 100, 5000, 5600),
    "Bp_m"
)

h_Bp_mass_kkpi.Draw("HIST")
c2.Update()
c2.Draw()


# In[9]:


# --- Salvar em arquivo ROOT ---
#out_file = ROOT.TFile("Bp_KKpi_hist.root", "RECREATE")
#h_Bp_mass_kkpi.Write("h_Bp_mass_kkpi")  # nome do histograma dentro do arquivo
#out_file.Close()


# In[10]:


from ROOT import (
    TFile, RooRealVar, RooDataHist, RooExponential, RooFit, TCanvas, kRed, RooArgList
)

# --- Abrir histograma ---
f = TFile.Open("Bp_KKpi_hist.root")
hist = f.Get("h_Bp_mass_kkpi")

# --- Variável de massa em MeV ---
x = RooRealVar("x", "Massa [MeV]", 5150, 5500)
x.setBins(hist.GetNbinsX())

# --- Dataset RooFit ---
data = RooDataHist("data", "Dados do histograma", RooArgList(x), RooFit.Import(hist))

# --- Modelo exponencial ---
tau = RooRealVar("tau", "tau", -0.001, -1.0, 0.0)
expo = RooExponential("expo", "Exponential PDF", x, tau)

# --- Definir ranges só para o fundo (fora do pico ~5280 MeV) ---
x.setRange("fundo_esq", 5196, 5255)   # região esquerda
x.setRange("fundo_dir", 5340, 5510)   # região direita

# --- Ajustar só ao fundo ---
expo.fitTo(data, RooFit.Range("fundo_esq,fundo_dir"))

# --- Plot em TODO o intervalo ---
c1 = TCanvas("c1", "Fit do fundo exponencial", 900, 700)
frame = x.frame(RooFit.Title("Fundo ajustado com exponencial"))
data.plotOn(frame)
expo.plotOn(frame, RooFit.LineColor(kRed), RooFit.Range("full"), RooFit.NormRange("fundo_esq,fundo_dir"))

frame.GetYaxis().SetTitle("Eventos")
frame.GetXaxis().SetTitle("M_{#pi^{+}K^{-}K^{+}} [MeV]")
frame.Draw()



c1.Update()
c1.Draw()


# In[ ]:





# In[ ]:





# In[ ]:




