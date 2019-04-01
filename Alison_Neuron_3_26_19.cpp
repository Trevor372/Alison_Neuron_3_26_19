#include <Alison_Neuron_3_26_19.h>
#include <math.h>

// Model Functions


//GAMMAF(VV,theta,sigma)=1.0/(1.0+exp(-(VV-theta)/sigma))
static inline double
m_inf(double V)
{
    //minf(V)=GAMMAF(V,thetam,sigmam)
    //1.0/(1.0+exp(-(V-thetam)/sigmam))
    return 1.0/(1.0+exp(-(V-(-36))/2));
}

static inline double
h_inf(double V)
{
    //hinf(V)=GAMMAF(V,thetah,sigmah)
    //1.0/(1.0+exp(-(V-thetah)/sigmah))
    return 1.0/(1.0+exp(-(V-(-45))/-3.0));
}

static inline double
n_inf(double V)
{
    //ninf(V)=GAMMAF(V,thetan,sigman)
    //1.0/(1.0+exp(-(V-thetan)/sigman))
    return 1.0/(1.0+exp(-(V - (-40.0))/5.0));
}

static inline double
tau_m(double V)
{
    //taum(V)=0.05+0.5*(GAMMAF(V,thetatma,sigmatma)*GAMMAF(V,thetatmb,sigmatmb))
    //taum(V)=0.05+0.5*((1.0/(1.0+exp(-(V-thetatma)/sigmatma))*1.0/(1.0+exp(-(V-thetatmb)/sigmatmb)))
    return 0.05+0.5*((1.0/(1.0+exp(-(V - (-20.0))/-10.0)))*(1.0/(1.0+exp(-(V - (-60.0))/3.0))));
}

static inline double
tau_h(double V, double tauh0, double tauh1)
{
    //tauh(V)=tauh0+(tauh1-tauh0)*GAMMAF(V,thetath,sigmath)
    //tauh(V)=tauh0+(tauh1-tauh0)*(1.0/(1.0+exp(-(V-thetath)/sigmath)))
    return tauh0+(tauh1-tauh0)*(1.0/(1.0+exp(-(V - (-45.0))/-3.0)));
}

static inline double
tau_n(double V, double taun0, double taun1)
{
    //taun(V)=taun0+(taun1-taun0)*GAMMAF(V,thetatna,sigmatna)
    //taun(V)=taun0+(taun1-taun0)*1.0/(1.0+exp(-(V-thetatna)/sigmatna))
    return taun0+(taun1-taun0)*(1.0/(1.0+exp(-(V - (-50.0))/-10.0)));
}

static inline double
a_inf(double V, double thetaa, double sigmaa)
{
    //ainf(v)=1/(1+exp(-(v-thetaa)/sigmaa))
    return 1.0/(1.0+exp(-(V - thetaa)/sigmaa));
}

static inline double
b_inf(double V, double thetab, double sigmab)
{
    //binf(v)=1/(1+exp(-(v-thetab)/sigmab))
    return 1.0/(1.0+exp(-(V - thetab)/sigmab));
}


extern "C" Plugin::Object *
createRTXIPlugin(void)
{
    return new Alison_Neuron_3_26_19();
}

static DefaultGUIModel::variable_t vars[] =
{
    { "Vm", "Membrane Potential", DefaultGUIModel::OUTPUT, },
    { "IA_INPUT", "Input from GA_Calc", DefaultGUIModel::INPUT, },
    { "IA-ACTIVATION", "Input from ia-activate", DefaultGUIModel::INPUT, },
    { "Iapp (nA)", "Applied Current",
        DefaultGUIModel::PARAMETER, },
    { "V0 (mV)", "Initial membrane potential (mV)",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "cm (nF)", "Specific membrane capacitance",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "IA_Scale", "Scalar multiplier of the IA current",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "Specific_GA", "Constant value determined from data",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "G_Na_max (uS)", "Maximum Na+ conductance density",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "E_Na (mV)", "Sodium reversal potential",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "G_K_max (uS)",
        "Maximum delayed rectifier conductance density",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },        
    { "E_K (mV)", "K+ reversal potential", DefaultGUIModel::PARAMETER
        | DefaultGUIModel::DOUBLE, },
    { "G_L (uS)", "Maximum leak conductance density",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "E_L (mV)", "Leak reversal potential",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "Rate (Hz)", "Rate of integration (Hz)", DefaultGUIModel::PARAMETER
        | DefaultGUIModel::UINTEGER, },
{ "thetaa", " ",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
{ "thetab", " ",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
{ "sigmaa", " ",DefaultGUIModel::PARAMETER | 
DefaultGUIModel::DOUBLE, },
{ "sigmab", " ",DefaultGUIModel::PARAMETER |
DefaultGUIModel::DOUBLE, },
    { "taua", "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "taub", "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	 { "taun0", "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "taun1", "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "tauh0", "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "tauh1", "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "m", "", DefaultGUIModel::STATE, },
    { "h", "", DefaultGUIModel::STATE, },
    { "n", "", DefaultGUIModel::STATE, },
    { "a", "", DefaultGUIModel::STATE, },
    { "b", "", DefaultGUIModel::STATE, },
    { "IA", "", DefaultGUIModel::STATE, },
    { "Time (s)", "Time (s)", DefaultGUIModel::STATE, }, };

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

/*
* Macros for making the code below a little bit cleaner.
*/

#define V (y[0])
#define m (y[1])
#define h (y[2])
#define n (y[3])
#define ns (y[4])
#define b (y[5])
#define a (y[6])
#define dV (dydt[0])
#define dm (dydt[1])
#define dh (dydt[2])
#define dn (dydt[3])
#define dns (dydt[4])
#define db (dydt[5])
#define da (dydt[6])
#define G_Na (G_Na_max*m*m*m*h)
#define G_K  (G_K_max*n*n*n*n)
#define G_Ks (G_Ks_max*ns*ns*ns*ns)
#define G_A  (Specific_GA*cm*IA_Scale*a*a*a*b)

Alison_Neuron_3_26_19::Alison_Neuron_3_26_19(void) : DefaultGUIModel("Alison_Neuron_3_26_19", ::vars, ::num_vars) {
    setWhatsThis(
    "<p><b>Alison_Neuron_3_26_19:</b><br>This module simulates a neuron with an A Current.</p>");
    createGUI(vars, num_vars);
    initParameters();
    update( INIT );
    refresh();
    resizeMe();
}

Alison_Neuron_3_26_19::~Alison_Neuron_3_26_19(void) {}

void Alison_Neuron_3_26_19::execute(void) {
    systime = count * period; // time in seconds
    for (int i = 0; i < steps; ++i){
        solve(period / steps, y); // period in s
	output(0) = V * 1e-3; // convert to V
    } 
    count++;
}

void Alison_Neuron_3_26_19::update(DefaultGUIModel::update_flags_t flag) {
    switch (flag) {
        case INIT:
            setParameter("V0 (mV)", QString::number(V0)); 
            setParameter("IA_Scale", QString::number(IA_Scale));
	    setParameter("Specific_GA", QString::number(Specific_GA));	
            setParameter("cm (nF)", QString::number(cm)); 
            setParameter("G_Na_max (uS)", QString::number(G_Na_max)); 
            setParameter("E_Na (mV)", QString::number(E_Na)); 
            setParameter("G_K_max (uS)", QString::number(G_K_max)); 
            setParameter("E_K (mV)", QString::number(E_K)); 
            setParameter("G_L (uS)", QString::number(G_L)); 
            setParameter("E_L (mV)", QString::number(E_L));  
            setParameter("Iapp (nA)", QString::number(Iapp)); 
            setParameter("Rate (Hz)", rate);
            setParameter("taua", QString::number(taua));
	    setParameter("taub", QString::number(taub));
	    setParameter("thetaa", QString::number(thetaa));
	    setParameter("thetab", QString::number(thetab));
            setParameter("sigmaa", QString::number(sigmaa));
	    setParameter("sigmab", QString::number(sigmab));
	    setParameter("taun0", QString::number(taun0));
	    setParameter("taun1", QString::number(taun1));
            setParameter("tauh0", QString::number(tauh0));
            setParameter("tauh1", QString::number(tauh1));
            setState("m", m);
            setState("h", h);
            setState("n", n);
            setState("b", b);
            setState("a", a);
            setState("IA",IA);
            setState("Time (s)", systime);
            break;

        case MODIFY:
            V0 = getParameter("V0 (mV)").toDouble();
            cm = getParameter("cm (nF)").toDouble();
	    IA_Scale = getParameter("IA_Scale").toDouble();
	    Specific_GA = getParameter("Specific_GA").toDouble();
            G_Na_max = getParameter("G_Na_max (uS)").toDouble();
            E_Na = getParameter("E_Na (mV)").toDouble();
            G_K_max = getParameter("G_K_max (uS)").toDouble();
            E_K = getParameter("E_K (mV)").toDouble();
            G_L = getParameter("G_L (uS)").toDouble();
            E_L = getParameter("E_L (mV)").toDouble();
            Iapp = getParameter("Iapp (nA)").toDouble(); 
            rate = getParameter("Rate (Hz)").toDouble();
            taua = getParameter("taua").toDouble();
            taub = getParameter("taub").toDouble();
	    thetaa = getParameter("thetaa").toDouble();              
            thetab = getParameter("thetab").toDouble();
	    sigmaa = getParameter("sigmaa").toDouble();
            sigmab = getParameter("sigmab").toDouble();
            taun0 = getParameter("taun0").toDouble();
            taun1 = getParameter("taun1").toDouble();
            tauh0 = getParameter("tauh0").toDouble();
            tauh1 = getParameter("tauh1").toDouble();
            steps = static_cast<int> (ceil(period * rate));

            break;

        case PERIOD:
            period = RT::System::getInstance()->getPeriod() * 1e-6; // time in seconds
            steps = static_cast<int> (ceil(period * rate));
            break;

        default:
            break;
    }
}

void Alison_Neuron_3_26_19::initParameters() {
    V0 = -55.038; // mV
    G_Na_max = 0.24;
    G_K_max = 0.011;
    G_L = 0.0018;
    E_Na = 50.0; // mV
    E_K = -100.0;
    E_L = -59.5;
    cm = 0.0187;
    IA_Scale = 1;
    Specific_GA = 0.438;
    Iapp = 0.0; // 1 Hz spiking
    rate = 400;
    taua = 2.0;
    taub = 140.0;
    thetaa = -52.3;
    thetab = -79.0;
    sigmaa = 22.3;
    sigmab = -3.0;	
    taun0 = 1.8;
    taun1 = 138;
    tauh0 = 1.5;
    tauh1 = 135;
    V = V0;
    m = 0.00001;
    h = 0.8522;
    n = 0.000208;
    a = a_inf(V0, thetaa, sigmaa);
    b = b_inf(V0, thetab, sigmab);
    count = 0;
    systime = 0;
    period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
    steps = static_cast<int> (ceil(period * rate)); // calculate how many integrations to perform per execution step
}

void Alison_Neuron_3_26_19::solve(double dt, double *y) {
    double dydt[7];
    derivs(y, dydt);
    for (size_t i = 0; i < 7; ++i){
        y[i] += dt * dydt[i];
    }
}

void Alison_Neuron_3_26_19::derivs(double *y, double *dydt) {

    if(IA_Scale == 0){
	IA = -(input(0)); 
    }
    else{
	IA = G_A * (V - E_K);
    }	

    dV = (-(G_Na*(V-E_Na) + G_K*(V-E_K) + G_L*(V-E_L) + IA) + (Iapp + input(1)))/cm;
    dm = (m_inf(V) - m) / tau_m(V);
    dh = (h_inf(V) - h) / tau_h(V,tauh0,tauh1);
    dn = (n_inf(V) - n) / tau_n(V,taun0,taun1);
    db = (b_inf(V, thetab, sigmab) - b) / taub;
    da = (a_inf(V, thetaa, sigmaa) - a) / taua;
}
