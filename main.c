// MIT License

// Copyright (c) 2025 Martin Novak

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <stdio.h>
#include <math.h>
#include <windows.h>
#include <time.h>

#define resolucijaX 1900
#define resolucijaY 650

#define sirinaHeightmap 1000
#define dolzinaHeightmap 1000

// Barva ozadja (iz wincon - windows console knjižnice)
#define barvaOzadja (WORD)0xFFFF

// Ravnina na katero "slikam".
#define projekcijskaRavnina 120.0f

// Ravnina na katero se clippa trikotnike.
#define nearClip 0.1f

// Ravnina, ki odreže, kar je vidno a predaleč od "kamere" (nujno za normalen performance pri velikih heightmapih).
#define farClip 2000.0f

#define amplitudaTerena 300

// Koeficient razmaka med vsakim indeksom v heightmap arraju.
#define koefRazdalje 30.0f

// Buffer ki se bo izpisal direktno v terminal.
CHAR_INFO buffer[resolucijaY][resolucijaX];

// Buffer preko katerega bom ločeval med ploskvami (edino tiste najbližje "kameri" se bodo izrisale).
float ZBuffer[resolucijaY][resolucijaX]={0.0f};

// Array z informacijami za heightmap (višino pri posamezni poziciji).
float heightmap[dolzinaHeightmap][sirinaHeightmap];

// Max število barv je 16, tako da sem tu dokaj limitiran. Uporabljam pa barve od
// modre do rdeče, kjer so modri trikotniki "najnižji".
WORD barve[6]={
    BACKGROUND_RED | BACKGROUND_INTENSITY,
    BACKGROUND_RED | BACKGROUND_GREEN, 
    BACKGROUND_RED | BACKGROUND_GREEN | BACKGROUND_INTENSITY,
    BACKGROUND_GREEN | BACKGROUND_INTENSITY,  
    BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_INTENSITY,  
    BACKGROUND_BLUE                                                                           
};

// Struct za točke.
typedef struct {
    float x, y, z;
} vektor3D;

// Struct za trikotnik.
typedef struct {
    vektor3D oglisce[3];
    vektor3D normala;
} trikotnik;

// Struct za hranjenje informacij o rotaciji.
typedef struct{
    float kotX;
    float kotY;
    float kotZ;
} rotacija;

// Precomputane vrednosti sin/cos.
typedef struct{
    float sinX, cosX;
    float sinY, cosY;
    float sinZ, cosZ;
} precompTrig;

// Transformacije posameznega objekta/globalne transformacije.
typedef struct{
    rotacija ROT;
    float razmerje;
    vektor3D premik;
} transformacije;

// Struct za ploskve modela.
typedef struct{
    size_t velikost;
    size_t idxCURR;                 // Trenutni indeks.
    transformacije lastnosti;
    trikotnik ogliscaARR[];
} model;

// Inicializiramo globalne transformacije.
transformacije globalneTransformacije;

// Za debugganje vektorjev.
void printfV(vektor3D t){
    printf("%f, %f, %f\n", t.x, t.y, t.z);
}

void pocistiModel(model* modelStruct){modelStruct->idxCURR=0;}

void generirajTeren(){

    // Nastavim seed za RNG.
    srand(time(NULL));

    // Nastavim višine za vsak indeks v HM arraju. Delim z amplitudo terena, da dobim omejene višine.
    for (int x=0; x<dolzinaHeightmap; x++){
        for (int z=0; z<sirinaHeightmap; z++){
            heightmap[x][z]=(float)(abs(rand()%amplitudaTerena));
        }
    }

    // V tem delu sem se zgledoval po CNN nevronskih mrežah, kjer se uporablja "okno", ki vrne povprečje vseh sosednjih pikslov.
    for (int x=1; x<dolzinaHeightmap-1; x++){
        for (int z=1; z<sirinaHeightmap-1; z++){

            // Vse sosednje točke seštejemo in na koncu delimo z 9.0f, da dobimo povprečje.
            float vsota=0.0f;
            for (int dx=-1; dx<=1; dx++){
                for (int dz=-1; dz<=1; dz++){
                    vsota+=heightmap[x+dx][z+dz];
                }
            }

            heightmap[x][z]=vsota/9.0f;
        }
    }
}

// Inicializira model in alocira hardlimitano število ploskev.
model* initModel(size_t velikostARR){

    // Alociram spomin za osnovno strukturo model, in dodatnih 'velikostARR' elementov tipa trikotnik
    model* modelStruct=malloc(sizeof(model)+velikostARR*sizeof(trikotnik));

    // ERRORNE in exitne.
    if (!modelStruct){
        printf("Ni bilo mogoce mallocati.\n");
        exit(1);
    }

    modelStruct->velikost=velikostARR;
    modelStruct->idxCURR=0;

    return modelStruct;
}

// Push-ne ploskev na array z ploskvami.
int pushTriktonik(model* modelStruct,  trikotnik el){

    if (modelStruct->idxCURR==modelStruct->velikost) return -1;

    modelStruct->ogliscaARR[modelStruct->idxCURR]=el;
    modelStruct->idxCURR++;;

    return 0;
}

// Vse Vektorske funkcije se samo pokličejo na začetku, ko se vsi trikotniki naložijo, da avtomatsko precomputa vse normale itd... Lahko pride prav pri kakih svetlobnih izračunih.
//====================================================================================================================================================================================
inline float dolzina(vektor3D a){
    return sqrtf(a.x*a.x+a.y*a.y+a.z*a.z);
}

vektor3D normaliziraj(vektor3D a){
    
    float d=dolzina(a);
    
    if (d==0) return (vektor3D){0.0f, 0.0f, 0.0f};
    
    return (vektor3D){a.x/d, a.y/d, a.z/d};
}

void vektorskiProdukt(const vektor3D a, const vektor3D b, const vektor3D c, vektor3D* normala){

    vektor3D u={b.x-a.x, b.y-a.y, b.z-a.z};
    vektor3D v={c.x-a.x, c.y-a.y, c.z-a.z};

    *normala=(vektor3D){u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x};
    *normala=normaliziraj(*normala);
}

inline float skalarniProdukt(const vektor3D a, const vektor3D b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}  
// Vse Vektorske funkcije se samo pokličejo na začetku, ko se vsi trikotniki naložijo, da avtomatsko precomputa vse normale itd... Lahko pride prav pri kakih svetlobnih izračunih.
//====================================================================================================================================================================================

void nastaviTrikotnik(const vektor3D a, const vektor3D b, const vektor3D c, trikotnik *trikotnik){
    
    trikotnik->oglisce[0]=a;
    trikotnik->oglisce[1]=b;
    trikotnik->oglisce[2]=c;

    vektorskiProdukt(a, b, c, &trikotnik->normala);
}

void razcleniHeightmap(model* objekt){

    trikotnik TMP;

    // Offset, da lahko postavim mapo v sredino kooordinatnega sistema. Pomnoženo s koeficientom, zaradi razdalje med posameznim indeksom v tabeli za heightmap.
    float offx=(dolzinaHeightmap/2)*koefRazdalje;
    float offy=(sirinaHeightmap/2)*koefRazdalje;

    // Trikotnike se mora v pravilnem vrstnem redu naložiti, sicer se normale ne izračunajo pravilno.
    for (size_t idx=0; idx<dolzinaHeightmap-1; idx++){
        for (size_t idy=0; idy<sirinaHeightmap-1; idy++){
            
            nastaviTrikotnik(   (vektor3D){koefRazdalje*idx-offx,
                                        heightmap[idx][idy],
                                        koefRazdalje*idy-offy},

                                (vektor3D){koefRazdalje*idx-offx,
                                        heightmap[idx][idy+1],
                                        koefRazdalje*(idy+1)-offy}, 

                                (vektor3D){koefRazdalje*(idx+1)-offx,
                                        heightmap[idx+1][idy],
                                        koefRazdalje*idy-offy},
                                &TMP);

            pushTriktonik(objekt, TMP);


            nastaviTrikotnik(   (vektor3D){koefRazdalje*(idx+1)-offx,
                                        heightmap[idx+1][idy+1],
                                        koefRazdalje*(idy+1)-offy}, 

                                (vektor3D){koefRazdalje*(idx+1)-offx,
                                        heightmap[idx+1][idy],
                                        koefRazdalje*idy-offy}, 

                                (vektor3D){koefRazdalje*idx-offx,
                                        heightmap[idx][idy+1],
                                        koefRazdalje*(idy+1)-offy}, 
                                &TMP);

            pushTriktonik(objekt, TMP);
        }
    }

    objekt->lastnosti=(transformacije){{0, 0, 0}, 1.0f, {0, 0, 0}};
}

inline vektor3D premakni(vektor3D vektor, vektor3D premik){ 
    
    return (vektor3D){
        vektor.x+premik.x,
        vektor.y+premik.y,
        vektor.z+premik.z
    };
}

inline vektor3D rotiraj(vektor3D vektor, const precompTrig* rotacije){
    
    return (vektor3D){
        rotacije->cosY*rotacije->cosZ*vektor.x+(rotacije->cosY*rotacije->sinZ*vektor.y)-(rotacije->sinY*vektor.z), 
        (rotacije->sinX*rotacije->sinY*rotacije->cosZ-rotacije->cosX*rotacije->sinZ)*vektor.x+(rotacije->sinX*rotacije->sinY*rotacije->sinZ+rotacije->cosX*rotacije->cosZ)*vektor.y+rotacije->sinX*rotacije->cosY*vektor.z,
        (rotacije->cosX*rotacije->sinY*rotacije->cosZ+rotacije->sinX*rotacije->sinZ)*vektor.x+(rotacije->cosX*rotacije->sinY*rotacije->sinZ-rotacije->sinX*rotacije->cosZ)*vektor.y+rotacije->cosX*rotacije->cosY*vektor.z
    };
}

inline vektor3D perspektiva(vektor3D vektor){

    return (vektor3D){
        (vektor.x*projekcijskaRavnina)/(vektor.z),
        (vektor.y*projekcijskaRavnina)/(vektor.z),
        (vektor.z)
    };
}

void transformacijaTrikotnika(trikotnik* tempTrikotnik, const trikotnik* originalniTrikotnik, const precompTrig* lokalneRotacije){
    
    tempTrikotnik->normala=rotiraj(originalniTrikotnik->normala, lokalneRotacije);

    for (size_t idx=0; idx<3; idx++) tempTrikotnik->oglisce[idx]=premakni(rotiraj(originalniTrikotnik->oglisce[idx], lokalneRotacije), globalneTransformacije.premik);
}

// Iterira skozi celotni buffer in ga "počisti".
void pocistiBuffer(){
    
    for (size_t y=0; y<resolucijaY; y++){
        for (size_t x=0; x<resolucijaX; x++){
            buffer[y][x].Attributes=barvaOzadja;
            ZBuffer[y][x]=INFINITY;
        }
    }
}

// Postavi piksle, obene pregleda če je točka znotraj zaslona ali ni, potem še preveri z-buffer da naprej nariše bližnje ploskve.
void postaviPixel(int x, int y, float z, WORD barva){
    
    if (x+resolucijaX/2>=0 && x+resolucijaX/2<resolucijaX && y+resolucijaY/2>=0 && y+resolucijaY/2<resolucijaY && (z<ZBuffer[y+resolucijaY/2][x+resolucijaX/2])){
        buffer[y+resolucijaY/2][x+resolucijaX/2].Attributes=barva;
        ZBuffer[y+resolucijaY/2][x+resolucijaX/2]=z;
    }
}

// https://learn.microsoft.com/en-us/windows/console/writeconsoleoutput. Lažje manipulacije v terminalu brez ansi escape kod.
void renderajBuffer(){
    
    SMALL_RECT okno={0, 0, resolucijaX-1, resolucijaY-1};
    WriteConsoleOutput(GetStdHandle(STD_OUTPUT_HANDLE), (CHAR_INFO*)buffer, (COORD){resolucijaX, resolucijaY}, (COORD){0, 0}, &okno);
}

// Risanje trikotnika s težiščnimi koordinatami. Risanje teh ni trivialno, zato sem funkcijo (sicer modificirano) napisal po psevdokodi z
// https://www.sunshine2k.de/coding/java/TriangleRasterization/TriangleRasterization.html.
void zapolniTrikotnik(trikotnik *t, WORD barva){

    // Zamenja točke po višini.
    if (t->oglisce[1].y<t->oglisce[0].y) {vektor3D tmp=t->oglisce[0]; t->oglisce[0]=t->oglisce[1]; t->oglisce[1]=tmp;}
    if (t->oglisce[2].y<t->oglisce[0].y) {vektor3D tmp=t->oglisce[0]; t->oglisce[0]=t->oglisce[2]; t->oglisce[2]=tmp;}
    if (t->oglisce[2].y<t->oglisce[1].y) {vektor3D tmp=t->oglisce[1]; t->oglisce[1]=t->oglisce[2]; t->oglisce[2]=tmp;}

    // Ker uporabljam izredno nestandarden način renderanja (obenem se večina teh zadev ne dela več ročno pri dejanskih graphics api-jih),
    // je prišlo do bug-a, kjer se je celotni program zrušil ali je deloval izredno počasno (cca 5s na frame), kadar je bilo ogromno clippanja.
    // To je bilo zato, ker prejšnja funkcija ne upošteva, da pri projekciji na ravnino, če se točka nahaja za njo, vrednosti x in y bistveno hitreje
    // naraščajo, zato so te vrednosti "eksplodirale". S kombinacijami fmax(fmin(fmin())) in fmin(fmax(fmax())) te vrednosti omejim.
    int minX=(int)fmax(fminf(fminf(t->oglisce[0].x, t->oglisce[1].x), t->oglisce[2].x), -resolucijaX/2.0f);
    int maxX=(int)fminf(fmaxf(fmaxf(t->oglisce[0].x, t->oglisce[1].x), t->oglisce[2].x), resolucijaX/2.0f);
    int minY=(int)fmaxf(fminf(fminf(t->oglisce[0].y, t->oglisce[1].y), t->oglisce[2].y), -resolucijaY/2.0f);
    int maxY=(int)fminf(fmaxf(fmaxf(t->oglisce[0].y, t->oglisce[1].y), t->oglisce[2].y), resolucijaY/2.0f);

    vektor3D v0=t->oglisce[0];
    vektor3D v1=t->oglisce[1];
    vektor3D v2=t->oglisce[2];

    // To je 2x površine. (izracun preko determinante).
    float povrsina=(v1.y-v2.y)*(v0.x-v2.x)+(v2.x-v1.x)*(v0.y-v2.y);
    if (fabsf(povrsina)<1e-6f) return;

    // Pregledamo vse min/max vrednosti koordinat, da vidimo če se točke ujemajo s pogojem, kjer morajo biti l1, l2, l3 vsi >= 0.
    // Pogoj je v razdelku "Pretvorba v težiščne koordinate" na naslovu https://sl.wikipedia.org/wiki/Te%C5%BEi%C5%A1%C4%8Dni_koordinatni_sistem
    for (int y=minY; y<=maxY; y++){
        for (int x=minX; x<=maxX; x++){
            float px=(float)x;
            float py=(float)y;

            float l1=((v1.y-v2.y)*(px-v2.x)+(v2.x-v1.x)*(py-v2.y))/povrsina;
            float l2=((v2.y-v0.y)*(px-v2.x)+(v0.x-v2.x)*(py-v2.y))/povrsina;
            float l3=1.0f-l1-l2;

            if (l1>=0 && l2>=0 && l3>=0){
                float z=l1*v0.z+l2*v1.z+l3*v2.z;
                postaviPixel(x, y, z, barva);
            }
        }
    }
}

// Helper funkcija, ki išče presečišče stranice trikotnika s projekcijsko ravnjo.
inline vektor3D presecisce(vektor3D a, vektor3D b, float t){
    
    return (vektor3D){
        a.x+(b.x-a.x)*t,
        a.y+(b.y-a.y)*t,
        a.z+(b.z-a.z)*t
    };
}

void clippanje(trikotnik* ploskev, WORD barva){
    
    int stevecNotranje=0, stevecZunanje=0;
    int notranjeIDX[3], zunanjeIDX[3];
    static trikotnik tempTri1, tempTri2;
    static vektor3D t1, t2;

    // Razvrsti oglišča glede na nearClip.
    for (int i=0; i<3; i++){
        if (ploskev->oglisce[i].z>=nearClip) notranjeIDX[stevecNotranje++]=i;
        else zunanjeIDX[stevecZunanje++]=i;
    }

    // Vse zunaj.
    if (stevecNotranje==0) return;

    // Ena znotraj nastane en trikotnik
    else if (stevecNotranje==1){
        vektor3D n=ploskev->oglisce[notranjeIDX[0]];
        vektor3D z1=ploskev->oglisce[zunanjeIDX[0]];
        vektor3D z2=ploskev->oglisce[zunanjeIDX[1]];

        float k1=(nearClip-n.z)/(z1.z-n.z);
        float k2=(nearClip-n.z)/(z2.z-n.z);

        t1=presecisce(n, z1, k1);
        t2=presecisce(n, z2, k2);

        tempTri1.oglisce[0]=perspektiva(n);
        tempTri1.oglisce[1]=perspektiva(t1);
        tempTri1.oglisce[2]=perspektiva(t2);

        zapolniTrikotnik(&tempTri1, barva);
    }

    // Dve točki znotraj, torej nastaneta dva trikotnika
    else if (stevecNotranje==2){
        vektor3D n1=ploskev->oglisce[notranjeIDX[0]];
        vektor3D n2=ploskev->oglisce[notranjeIDX[1]];
        vektor3D z=ploskev->oglisce[zunanjeIDX[0]];

        float k1=(nearClip-n1.z)/(z.z-n1.z);
        float k2=(nearClip-n2.z)/(z.z-n2.z);

        t1=presecisce(n1, z, k1);
        t2=presecisce(n2, z, k2);

        
        tempTri1.oglisce[0]=perspektiva(n1);
        tempTri1.oglisce[1]=perspektiva(n2);
        tempTri1.oglisce[2]=perspektiva(t1);

        tempTri2.oglisce[0]=perspektiva(n2);
        tempTri2.oglisce[1]=perspektiva(t1);
        tempTri2.oglisce[2]=perspektiva(t2);

        zapolniTrikotnik(&tempTri1, barva);
        zapolniTrikotnik(&tempTri2, barva);
    }

    // Vse tri znotraj.
    else{
        zapolniTrikotnik(&(trikotnik){perspektiva(ploskev->oglisce[0]), perspektiva(ploskev->oglisce[1]), perspektiva(ploskev->oglisce[2])}, barva);
    }
}

// Helper funkcija, ki precomputa sinuse in kosinuse.
void lokalniPrecompute(model* objekt, precompTrig* rotacije){
    
    rotacije->sinX=sinf(objekt->lastnosti.ROT.kotX);
    rotacije->cosX=cosf(objekt->lastnosti.ROT.kotX);
    rotacije->sinY=sinf(objekt->lastnosti.ROT.kotY);
    rotacije->cosY=cosf(objekt->lastnosti.ROT.kotY);
    rotacije->sinZ=sinf(objekt->lastnosti.ROT.kotZ);
    rotacije->cosZ=cosf(objekt->lastnosti.ROT.kotZ);
}

void rasterizacija(model* objekt){
    
    size_t d=objekt->velikost;

    // Vnaprej zračuna vse trig funkcije, da hitreje izvede rotacije.
    static precompTrig lokalneRotacije;
    static trikotnik t;

    lokalniPrecompute(objekt, &lokalneRotacije);

    WORD barva;

    for (size_t idx=0; idx<d; idx++){

        transformacijaTrikotnika(&t, &objekt->ogliscaARR[idx], &lokalneRotacije);
        
        // Backface culla (če je skalarni produkt < 0 se trikotnik ne vidi).
        float skalar=skalarniProdukt(t.normala, normaliziraj((vektor3D){t.oglisce[0].x, t.oglisce[0].y, t.oglisce[0].z}));
        if (skalar<0 || dolzina((vektor3D){t.oglisce[0].x, t.oglisce[0].y, t.oglisce[0].z})>farClip) continue;

        // 1. aritmetična sredina y koordinate -> delim z amplitudo s številom barv...
        else barva=barve[(int)(6*(objekt->ogliscaARR[idx].oglisce[0].y+objekt->ogliscaARR[idx].oglisce[1].y+objekt->ogliscaARR[idx].oglisce[2].y)/(3.0f*amplitudaTerena))];
        
        clippanje(&t, barva);
    }
    
    
}

void preverjanjeInputa(model* objekt){

    if (GetAsyncKeyState('W') & 0x8000){
        globalneTransformacije.premik.z-=2.5;
    }
    if (GetAsyncKeyState('S') & 0x8000){
        globalneTransformacije.premik.z+=2.5;
    }
    if (GetAsyncKeyState('F') & 0x8000){
        objekt->lastnosti.ROT.kotY+=0.01f;
    }
    if (GetAsyncKeyState('G') & 0x8000){
        objekt->lastnosti.ROT.kotY-=0.01f;
    }
    if (GetAsyncKeyState('E') & 0x8000){
        globalneTransformacije.premik.y+=3.0f;
    }
    if (GetAsyncKeyState('Q') & 0x8000){
        globalneTransformacije.premik.y-=3.0f;
    }
    if (GetAsyncKeyState('R') & 0x8000){
        pocistiModel(objekt);
        generirajTeren();
        razcleniHeightmap(objekt);
    }
    if (GetAsyncKeyState(VK_ESCAPE) & 0x8000){
        free(objekt);
        exit(0);
    }
}

int main(){

    globalneTransformacije=(transformacije){{0, 0, 0}, 1.0f, {0, 0, 450}};
    
    model* objekt=initModel((sirinaHeightmap-1)*(dolzinaHeightmap-1)*2); 

    generirajTeren();

    razcleniHeightmap(objekt);

    while (1){
        preverjanjeInputa(objekt);
        pocistiBuffer();
        rasterizacija(objekt);
        renderajBuffer();
    }

    return 0;

}
