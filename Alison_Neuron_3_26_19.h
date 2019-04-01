#include <default_gui_model.h>


class Alison_Neuron_3_26_19 : public DefaultGUIModel
{
    
    public:
    
        Alison_Neuron_3_26_19(void);
        virtual~Alison_Neuron_3_26_19(void);
    
        void execute(void);
    
    protected:
    
        void update(DefaultGUIModel::update_flags_t);
    
    private:
    
        void derivs(double *, double *);
        void solve(double, double *);
        void initParameters();
    
        double y[7];
        double period;
        int steps;
    
        double V0;
	double IA_Scale;
	double Specific_GA;
        double cm;
        double G_Na_max;
        double E_Na;
        double G_K_max;
        double E_K;
        double G_L;
        double E_L;
        double G_A_max;
        double Iapp;
        double rate;
        double taua;
	double taub;
	double thetaa;
	double thetab;
	double sigmaa;
	double sigmab;
	double taun0;
	double taun1;
	double tauh0;
	double tauh1;
        double IA;
        double systime;
        long long count;
};
