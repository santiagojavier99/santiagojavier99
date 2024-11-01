#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <complex>
#include <string>

using namespace std;

int main () {

int Nx=1000;
int Nts=1500;
double dt=0.1; 
const double Vo=4; 
const double a=250; 
const double b=260; 
double D=500;
double t=0;
int Nprint=1;
double normaq, norma;
complex<double> imm(0.0,1.0);


complex<double> *psi = new complex<double>[Nx];


double h=D/(Nx-1);
normaq=0.;
/*cout << "Inserire punti lungo x" << endl;
cin >> Nx;
cout << "Inserire condizione al contorno [psi(x)=psi(x+D)]" << endl;
cin >> D;
cout << "Inserire Vo" << endl;
cin >> Vo;
cout << "Inserire limite inferiore (a) e superiore (b) di barriera di potenziale" << endl;
cin >> a >> b;
cout << "Inserire dt" << endl;
cin >> dt;
cout << "Inserire numero di time steps " << endl;
cin >> Nts;
cout << "Scrivi ogni numero passi\n";
cin >> Nprint;
*/


char scelta;
cout << "Quale funzione si vuole avere?" << endl;
cin >> scelta;
switch(scelta)
{
case 'a': //seno
cout << "Hai scelto la funzione seno" << endl;
for(int i=1; i<Nx-1; i++)
	{
	double x=h*i;
	psi[i]=sin(2*M_PI*x/D); // Il 2*pi dà condizioni periodiche alla funzione d'onda
	normaq = normaq + pow(abs(psi[i]),2);
	}
break;

case 'b': //gaussiana
cout << "Hai scelto la funzione gaussiana" << endl;
double x_o=150;
double k_o=3;
long int N_int=2000;
double sigma=0.1; // Scelta di sigma dovuta a una dispersione del pacchetto.
const double hk=sigma*10/N_int; // Per controllare la precisione con cui faccio i passi-momento a seconda della sigma immessa

/*cout << "Per quanto riguarda il pacchetto gaussiano: " << endl;
cout << "Il pacchetto è centrato nella posizione " << endl;
cin >> x_o;
cout << "Il pacchetto è centrato nel momento " << endl;
cin >> k_o;
cout << "La sigma della gaussiana è " << endl;
cin >> sigma;*/

//Definisco array necessari per trattare il pacchetto gaussiano in trasformata
complex<double> *griglia = new complex<double>[2*N_int];
complex<double> *f = new complex<double>[2*N_int];
complex<double> *integrale = new complex<double>[2*N_int];
complex<double> somma1;
complex<double> *somma = new complex<double>[Nx];
double modulo=0.;

for(int i=1; i<Nx-1; i++)
	{
	double x=D/(Nx-1)*i; 
	for(int k=-N_int; k<N_int; k++)
		{
		griglia[k]=k_o+k*hk;
		//cout << k << " " << griglia[k] << endl;
		f[k]=1./(sqrt(2*M_PI))*exp(imm*(griglia[k])*(x-x_o))*exp(-pow((griglia[k]-k_o)/sigma,2));
		//cout << f[k] << endl;
		//file << griglia[k] << " " << real(f[k]) << endl;
		}
	
	somma1=0.;

	for(int k=-N_int+1; k<N_int-1; k+=2)
		{ 
		integrale[k]=(f[k+1]/3.+4.*f[k]/3.+f[k-1]/3.)*hk;
		somma1 = somma1 + integrale[k];
		}
		//modulo = modulo + pow(abs(somma),2);
		
		somma[i]=somma1;
		//file << x << " " << real(somma[i]) << endl;
	}
		//file.close();
		
//HO DEFINITO PACCHETTO DI GAUSS, È IN ARRAY SOMMA.




for(int i=1; i<Nx-1; i++)
	{
	double x=h*i;
	psi[i]=somma[i];
	normaq = normaq + pow(abs(psi[i]),2);
	}

break;}

// EVOLUZIONE TEMPORALE:

//Definisco gli array complessi necessari per l'evoluzione temporale
complex<double> *psi_norm = new complex<double>[Nx];
complex<double> *m11 = new complex<double>[Nx];
complex<double> *alpha = new complex<double>[Nx];
complex<double> *beta = new complex<double>[Nx];
complex<double> *F = new complex<double>[Nx];
complex<double> *ai = new complex<double>[Nx];

psi[0]=0;
psi[Nx-1]=0;	
normaq=normaq*D/(Nx-1);
norma = sqrt(normaq);
//cout << "La norma è " << norma << endl;

for(int i=0; i<Nx; i++) // Normalizzo la funzione d'onda.
	{
	psi_norm[i]=psi[i]/norma;
	psi[i]=psi_norm[i];
	}
	
double V;
double x;
for(int j=0; j<Nx; j++)
	{
	x=h*j;
	if((x<a) or (x>b)) V=0;
	else V=Vo;
	m11[j]=imm*4.*pow(h,2)/dt-2.*pow(h,2)*V-2.; 
	//cout << m11[j] << " " << V << endl;
	}
// Ho finito di calcolare la matrice M_ij

ofstream fileg;
string nomeg;
nomeg = (string) "psi_completa.dat";
fileg.open(nomeg);
fileg.precision(10);

for(int k=0; k<Nts; k++) // per fare progredire i time steps
	{
	
	for(int i=1; i<Nx-1; i++)
		{
		double x=h*i;
		double V;
		if(x<a or x>b) V=0;
		else V=Vo;
		F[i]=-psi[i+1]+2.*psi[i]-psi[i-1]+(imm*4.*pow(h,2)/dt)*psi[i]+2.*pow(h,2)*V*psi[i];
		//cout<<F[i]<<endl;
		F[0]=-psi[1]+2.*psi[0]+imm*4.*pow(h,2)/dt*psi[0]+2.*pow(h,2)*V*psi[0];
		F[Nx-1]=2.*psi[Nx-1]+psi[Nx-2]+(imm*4.*pow(h,2)/dt)*psi[Nx-1]+2.*pow(h,2)*V*psi[Nx-1];
		}

	alpha[0]=-m11[0];
	beta[0]=F[0];
	
	for(int w=1; w<Nx-1; w++)
		{
		alpha[w]=-pow(alpha[w-1],-1)-m11[w];
		beta[w]=F[w]+(beta[w-1])/(alpha[w-1]);
		//cout<<w<<"	"<<m11[w]<<"	"<<alpha[w]<<"	"<<beta[w]<<endl;
		}
	ai[Nx-1]=pow((1./alpha[Nx-2]+m11[Nx-1]),-1)*(F[Nx-1]+beta[Nx-2]/alpha[Nx-2]);
        
	for(int l=Nx-1; l>0; l--)
		{
		ai[l-1]=(ai[l]-beta[l-1])/(alpha[l-1]);
		psi[l-1]=ai[l-1];
		}
		
	if(k%Nprint==0)
		{
     		ofstream file;
        	string nome;
        	nome = (string) "psi"+to_string(k) +(string)".dat";
            	file.open(nome);
            	file.precision(10);
            	for(int i=0;i<Nx-1;i++)
            		{
               	double x=h*i;
                	file << x << "  " << pow(abs(psi[i]),2) << endl; // plotta real
                	fileg << x << "  " << pow(abs(psi[i]),2) << endl;
               
            		}

            	file.close();
            	fileg << '\n';

        	}    
	norma=0.;
	
	// Controllo se la norma è a 1 ad ogni passo di time step.
 	for(int i=1;i<Nx-1;i++)
 		{
            	norma= norma + pow(abs(psi[i]),2);
        	}
        norma=norma*D/(Nx-1);
        cout << "Passo " << k << " Norma " << norma << endl;
}

fileg.close();
free(psi);
free(F);
free(m11);
free(alpha);
free(beta);

return 0;}
