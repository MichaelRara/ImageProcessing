unit Unit1;

interface

uses
  ShellAPI,Generics.Collections ,Jpeg, Windows, Messages, SysUtils, Variants, Classes, Graphics,
  Controls, Forms, Dialogs, StdCtrls, FileCtrl, OpenTIFF, strUtils,
  ExtCtrls,Math, Grids, ComCtrls;

type
  TMainForm = class(TForm)
    DriveComboBox1: TDriveComboBox;
    DirectoryListBox1: TDirectoryListBox;
    FileListBox1: TFileListBox;
    TvorbaSite: TButton;
    Korelace: TButton;
    StringGrid1: TStringGrid;
    Interpolace: TButton;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Krok: TEdit;
    Polomer: TEdit;
    PrahHod: TEdit;
    Label7: TLabel;
    Label8: TLabel;
    Label9: TLabel;
    Label10: TLabel;
    KrokPomSit: TEdit;
    PolomerPomSit: TEdit;
    Label11: TLabel;
    Label12: TLabel;
    OkoliKorelace: TEdit;
    PolomerKor: TEdit;
    Label13: TLabel;
    BilinearniInterpolace: TButton;
    Posuvy: TButton;
    MeritkoMapyPosuvu: TButton;
    UlozPosuvy: TButton;
    MaxPol: TEdit;
    Label14: TLabel;
    Resetuj: TButton;
    ProgressBar1: TProgressBar;
    procedure TvorbaSiteClick(Sender: TObject);
    procedure KorelaceClick(Sender: TObject);
    procedure InterpolaceClick(Sender: TObject);
    procedure BilinearniInterpolaceClick(Sender: TObject);
    procedure PosuvyClick(Sender: TObject);
    procedure MeritkoMapyPosuvuClick(Sender: TObject);
    procedure UlozPosuvyClick(Sender: TObject);
    procedure KrokChange(Sender: TObject);
    procedure PolomerChange(Sender: TObject);
    procedure PrahHodChange(Sender: TObject);
    procedure KrokPomSitChange(Sender: TObject);
    procedure PolomerPomSitChange(Sender: TObject);
    procedure MaxPolChange(Sender: TObject);
    procedure PolomerKorChange(Sender: TObject);
    procedure OkoliKorelaceChange(Sender: TObject);
    procedure ResetujClick(Sender: TObject);

  private
    { Private declarations }
  public
    { Public declarations }
  end;

type
    TPoleWord=array[0..10000*10000] of Word;
    TPole2Word=array[0..10000*10000] of Word;
    PPoleWord=^TPoleWord;
    PPole2Word=^TPole2Word;
    PRGBTripleArray = ^TRGBTripleArray;
    TRGBTripleArray = array[0..10000*10000] of TRGBTriple;
    TSingleArray = Array[0..10000*10000] of single;
    PSingleArray=^TSingleArray;
    TIntegerArray = Array[0..10000*10000] of Integer;
    PIntegerArray=^TIntegerArray;
    GPPointsRecordKor=^GPointsRecordKor;
    GPointsRecordKor = record
                          I,J,Position:Integer;
                          PosunI, PosunJ:Single;
                       end;
    GPFVertex=^FakeVertex1;
    FakeVertex1 = record I,J:Integer;  end;
    GPColors=^GColours;
    GColours = record
                R,G,B:Byte
              end;
var
  MainForm: TMainForm;
  Indexes: TIndexes;
  Image16Array: PImage16Array;
  TIFFinfo: TTIFFinfo;
  Images: Array[0..1000] of PImage16Array;
  GPointsPosition:PIntegerArray;  //Pole bodù, které obsahuje pozice bodù, které leí v síti
  GPointsKor,GLPixelData:Array of TList;
  GU,GSize,GHeader,Gn,Gwidth,Gheight,GPrahHodnota,GPolomer,GKrok:Integer;
  GB2: array of PPoleWord;    //Obsahuje naètená data (obrazy)
  GPruObraz: PPoleWord;      //Obsahuje zprùmìrovanı obraz z dat
  //GDicVertex:TDictionary<integer,GPPointsRecordKor>; // slovník obsahující vrcholy interpolaèní sítì
  AGDicVertex: array of TDictionary<integer,GPPointsRecordKor>; // slovník obsahující vrcholy interpolaèní sítì
  GDicPosition:TDictionary<Integer,integer>;
  GDicBrightness:TDictionary<Integer,single>;
  GArrayOfNames:array of string;//pole na názvy vstupních obrázkù
  GNameOfFolder:string;
implementation

{$R *.dfm}
function Orientation(x1, y1, x2, y2, Px, Py: integer): ShortInt;
var
  Orin: integer;
begin
  Orin := (x2 - x1) * (Py - y1) - (Px - x1) * (y2 - y1);
  if Orin > 0 then Result := 1    (* Orientaion is to the right-hand side  *)
  else if Orin < 0 then Result := -1  (* Orientaion is to the left-hand side   *)
  else
    Result := 0;                  (* Orientaion is neutral if result is 0  *)
end;

function CheckEdge(x1,x2,y1,y2,I,J:integer):boolean;   //zkontroluje, zda bod leí na hranici trojúhelníku
var a,b,c,d:Integer;
begin
  if x1>=x2 then//uspoøádáme I-té sloky podle velikosti
    begin
      a:=x1; b:=x2;
    end
    else
    begin
      a:=x2; b:=x1;
    end;
    if y1>=y2 then//uspoøádáme J-té sloky podle velikosti
    begin
      c:= y1; d:=y2;
    end
    else
    begin
      c:=y2; d:=y1;
    end;
  if ((b<=I) and (I<=a)) and ((d<=J) and (J<=c)) then result:=true//bod leí na pøímce mezi uzly trojúhelnkové sítì
  else result:=false;
end;

function Det2(J,I,x1,y1,x2,y2:Integer):boolean; //pro kontrolu, jestli bod leí na pøímce dvou bodù
var a:integer;
begin
  a:=(J-x1) * (I-y2) - (I-y1) * (J-x2);
  if (a=0) and (CheckEdge(y1,y2,x1,x2,I,J)) then result:= true//pokud a = 0, pak tyto tøi body leí na stejné pøímce, pokud CheckEdge=true, pak bod J,I leí na úseèce mezi zbylími dvìma body
  else result:=false;
end;

function CheckPixel(I,J:integer; T:array of GPPointsRecordKor ):boolean;  //Ovìøí, jesti bod leí v trojúhelníku nebo na jeho hranici, pokud ano, vrátí true
var  Or1, Or2, Or3: ShortInt;
begin
  Or1    := Orientation(T[0]^.J, T[0]^.I, T[1]^.J, T[1]^.I, J,I);
  if Or1=0 then result:=CheckEdge(T[0]^.I,T[1]^.I, T[0]^.J,T[1]^.J,I,J)
  else
  begin
    Or2    := Orientation(T[1]^.J, T[1]^.I, T[2]^.J, T[2]^.I, J,I);
    if Or2=0 then result:=CheckEdge(T[1]^.I,T[2]^.I,T[1]^.J,T[2]^.J,I,J)
    else
    begin
      Or3    := Orientation(T[2]^.J, T[2]^.I, T[0]^.J, T[0]^.I, J,I);
      if Or3=0 then result:=CheckEdge(T[2]^.I,T[0]^.I,T[2]^.J,T[0]^.J,I,J)
      else Result := (Or1 = Or2) and (Or2 = Or3);
    end;
  end;
end;

function CheckFakeVertex(LFakeVertex:TList;I, J, R2,S2:integer):boolean;  //vrátí true, pokud jsme nalezli pixel nevhodnı k interpolaci, jinak vrátí false
var C:integer; F:boolean; FVertex:GPFVertex;
begin
  F:=false;
  for C:=0 to (LFakeVertex.Count-1) do      //projede list vrcholù nevhodnıch k interpolaci
  begin
    FVertex:=LFakeVertex.Items[C];
    if (FVertex^.I=I+R2) and (FVertex^.J=J+S2) then //pokud zjistím, e aktuální pixel je vrchol nevhodnı k interpolaci, pak F := false a pøejde se na další záznam v listu vrcholù
    begin
      F:=true;
      break;
    end;
  end;
  result:=F;
end;

function Gradient(height, width, Size: integer; Data2:PPoleWord;Grad2:PSingleArray):single;
var I,J:Integer;
    R:Int64;
    max:Single;
    Grad:Array of Array of Single;
begin
  R:=Size div 2; //R je poèet pixelù v obraze, na jeden pixel pøipadají 2 bajty
  SetLength(Grad,3);   //Definujeme matici o 3 øádcích
  SetLength(Grad[0],R);//První øádek má R sloupcù, je to derivace podle x
  SetLength(Grad[1],R);//Druhı øádek má R sloupcù, je to derivace podle y
  SetLength(Grad[2],R);//Tøetí øádek má R sloupcù. je norma gradientu
  max:=0;
  for J:=0 to height-1 do
     begin
      for I:=0 to width-1 do
        begin
          //Derivace podle x
          if (I>0) and (I<(width-1)) then Grad[0,I+J*width]:=(Data2^[J*width+I-1]-Data2^[J*width+I+1])/2
          else
            begin
                if(I=0) then Grad[0,I+J*width]:=Data2^[J*width+I+1]-Data2^[J*width+I]
                else Grad[0,I+J*width]:=Data2^[J*width+I]-Data2^[J*width+I-1];
            end;
          //Derivace podle y
          if(J>0) and (J<height-1) then Grad[1,I+J*width]:=(Data2^[(J-1)*width+I]-Data2^[(J+1)*width+I])/2
          else
            begin
                if(J=0) then Grad[1,I+J*width]:=Data2^[(J+1)*width+I]-Data2^[J*width+I]
                else Grad[1,I+J*width]:=Data2^[J*width+I]-Data2^[(J-1)*width+I];
            end;
          //Vıpoèet Normy
           Grad[2,I+J*width]:=sqrt(Grad[0,I+J*width]*Grad[0,I+J*width] + Grad[1,I+J*width]*Grad[1,I+J*width] );
           if Grad[2,I+J*width]>max then max:=Grad[2,I+J*width];
           Grad2^[J*width+I]:=Grad[2,I+J*width];
        end;
     end;
    //SetLength(Grad,3);   //Definujeme matici o 3 øádcích
    SetLength(Grad[0],0);//První øádek má R sloupcù, je to derivace podle x
    SetLength(Grad[1],0);//Druhı øádek má R sloupcù, je to derivace podle y
    SetLength(Grad[2],0);//Tøetí øádek má R sloupcù. je norma gradientu
    SetLength(Grad,0);   //Definujeme matici o 3 øádcích
    result:=max;
end;

function CreateBmp(width, height,format:integer):TBitmap;
var Bitmap:TBitmap;
begin
  Bitmap:=TBitmap.Create;
  Bitmap.Width:=width;
  Bitmap.Height:=height;
  case format of
  24: Bitmap.PixelFormat:=pf24bit;
  end;
  result:=Bitmap;
end;

function HelpVertex(Uk,width,I2,J2,R2,S2,n:integer):boolean;
var L:integer;
    Clear:boolean;
    PRecordKor:GPPointsRecordKor;
begin
  Clear:=True;
  for L:=0 to (Uk-1) do
  begin
    PRecordKor:=GPointsKor[n].Items[L];
    if (I2+R2)*width+(J2+S2)=PRecordKor^.Position then
    begin
      clear:=false;
      break;
    end;
  end;
  result:=Clear;
end;

function Interpolation(I,J:integer; T: array of GPPointsRecordKor;S:boolean ):single; Overload;
var z0,z1,z2,v1,v2,v3,v4,v5,v6:single;
begin
  if S= true then
  begin
    z0:=T[0]^.PosunI;z1:=T[1]^.PosunI;z2:=T[2]^.PosunI;
  end
  else
  begin
    z0:=T[0]^.PosunJ;z1:=T[1]^.PosunJ;z2:=T[2]^.PosunJ;
  end;
  v1:=(T[2]^.I-T[0]^.I ) * (z1-z0) * (J-T[0]^.J);//hotovo
  v2:=(I-T[0]^.I) * (T[1]^.J-T[0]^.J) * (z2-z0);//hotovo
  v3:=(I-T[0]^.I) * (z1-z0) * (T[2]^.J-T[0]^.J);//hotovo
  v4:=(J-T[0]^.J) * (T[1]^.I-T[0]^.I) * (z2-z0);//hotovo
  v5:=(T[1]^.J-T[0]^.J) * (T[2]^.I-T[0]^.I);//hotovo
  v6:=(T[1]^.I-T[0]^.I) * (T[2]^.J-T[0]^.J);//hotovo
  if (v5-v6=0) then
  showmessage(floattostr(v5-v6));
  result := z0+((v1+v2-v3-v4)/(v5-v6));
end;

function BiInterpolation(J,I:single;width,n:integer):word;
var F1,F2,F3,F4:single;//dílèí interpolované hodnoty
    a,b,c,d:integer;// souøadnice
begin
  a:=Floor(J); b:=Ceil(J); c:=Floor(I); d:=Ceil(I);
  if (a=b) and (c=d) then result:=GB2[n]^[c*width+a]
  else if (a=b) and (c<>d) then result:=word(round(GB2[n]^[c*width+a] + ((GB2[n]^[d*width+a]-GB2[n]^[c*width+a]) / (d-c) ) * (I-c))) //svislá interpolace
  else if (a<>b) and (c=d) then result:=word(round(GB2[n]^[c*width+a] + ((GB2[n]^[c*width+b]-GB2[n]^[c*width+a]) / (b-a) ) * (J-a))) //vodorovná intepolace
  else
  begin
    F1:=((b-J) / (b-a)) * ((d-I) / (d-c)) * GB2[n]^[c*width+a];
    F2:=((J-a) / (b-a)) * ((d-I) / (d-c)) * GB2[n]^[c*width+b];
    F3:=((b-J) / (b-a)) * ((I-c) / (d-c)) * GB2[n]^[d*width+a];
    F4:=((J-a) / (b-a)) * ((I-c) / (d-c)) * GB2[n]^[d*width+b];
    result:=word(round(F1+F2+F3+F4));
  end;
end;

function Interpolation(J,I,x1,y1,x2,y2:Integer;z1,z2:single):single;Overload; //spoète interpolaci na pøímce
var t:single;
begin
  if (x2-x1)<>0 then
    begin
      t:=(J-x1)/(x2-x1);
      result:= z1 + (z2-z1)*t ;
    end
  else
    begin
      t:=(I-y1)/(y2-y1);
      result:= z1 + (z2-z1)*t ;
    end;
end;

procedure SavePixel(I,J,width,n:integer; Triangle:array of GPPointsRecordKor; PPixel: GPPointsRecordKor);
var smer:boolean;
begin
  PPixel^.I:=I; PPixel^.J:=J; PPixel^.Position:=I*width+J;
  Smer:=true; //pokud Smer =true, pak poèítej interpolaci pro posuv ve smìru I
  PPixel^.PosunI:=Interpolation(I,J,Triangle,Smer);
  Smer:=false; //pokud Smer =false, pak poèítej interpolaci pro posuv ve smìru J
  PPixel^.PosunJ:=Interpolation(I,J,Triangle,Smer);
  GLPixeldata[n].Add(PPixel);
end;

function SavePixelAndCheck(I,J,width,n:integer; Triangle:array of GPPointsRecordKor;PPixel:GPPointsRecordKor):boolean;
begin
  PPixel^.I:=I; PPixel^.J:=J; PPixel^.Position:=I*width+J;
  PPixel^.PosunI:=Interpolation(J,I,Triangle[0]^.J,Triangle[0]^.I,Triangle[1]^.J,Triangle[1]^.I,Triangle[0]^.PosunI,Triangle[1]^.PosunI);
  PPixel^.PosunJ:=Interpolation(J,I,Triangle[0]^.J,Triangle[0]^.I,Triangle[1]^.J,Triangle[1]^.I,Triangle[0]^.PosunJ,Triangle[1]^.PosunJ);
  GLPixeldata[n].Add(PPixel);
  result:=true;
end;

function AritmetickyPrumer(r,I,J,width:integer):extended;  //vrátí aritmetickı prùmer hodnot pixelù v okolí uzlu
var K,G1,G2,G3,G4,t:integer;
    Sumx,v:extended;
begin
  Sumx:=GPruObraz^[I*width+J]; //Zjistí jas uzlu sítì a pøidá jej do sumy, uzel leí ve zprùmìrovaném obraze
  t:=3;v:=(2*r+1)*(2*r+1);//v je poèet pixelù v okolí uzlu sítì
  for K:=1 to r do // Cyklus pro danı poèet krunic kolem uzlu, poèet krunic= polomìr nìjvìtší krunice tj. r
    begin
      if (K>1) then t:=t+2; //Nastaví novou délku sloupce, resp. øádku
        for G1:=0 to (t-1) do //Projede spodní øádek
          Sumx:=Sumx+GPruObraz^[(I+K)*width+(J-K+G1)];  //Pokud bod leí ve správném okolí pøiète se jeho hodnota jasu k sumì
        for G2:=0 to (t-1) do //Projede horní øádek
          Sumx:=Sumx+GPruObraz^[(I-K)*width+(J-K+G2)];  //Pokud bod leí ve správném okolí pøiète se jeho hodnota jasu k sumì
        for G3:=1 to (t-2) do //Projede levı sloupec
          Sumx:=Sumx+GPruObraz^[(I-K+G3)*width+(J-K)];  //Pokud bod leí ve správném okolí pøiète se jeho hodnota jasu k sumì
        for G4:=1 to (t-2) do //Projede pravı sloupec
          Sumx:=Sumx+GPruObraz^[(I-K+G4)*width+(J+K)];  //Pokud bod leí ve správném okolí pøiète se jeho hodnota jasu k sumì
    end;
  result:=Sumx/v;//Aritmetickı prùmìr pro PruObraz
end;

function SmerodatnaOdchylka(r,I,J,width:integer;Ax:extended):extended;
var K,G1,G2,G3,G4,t,v:integer;
    Sx:extended;
begin
  t:=3;v:=(2*r+1)*(2*r+1);//v je poèet pixelù v okolí uzlu sítì
  Sx:=sqr(GPruObraz^[I*width+J]-Ax);
  for K:=1 to r do // Cyklus pro danı poèet krunic kolem uzlu, poèet krunic= polomìr nìjvìtší krunice tj. r
    begin
      if (K>1) then t:=t+2; //Nastaví novou délku sloupce, resp. øádku
        for G1:=0 to (t-1) do //Projede spodní øádek
          Sx:=Sx+sqr(GPruObraz^[(I+K)*width+(J-K+G1)]-Ax);
        for G2:=0 to (t-1) do //Projede horní øádek
          Sx:=Sx+sqr(GPruObraz^[(I-K)*width+(J-K+G2)]-Ax);
        for G3:=1 to (t-2) do //Projede levı sloupec
          Sx:=Sx+sqr(GPruObraz^[(I-K+G3)*width+(J-K)]-Ax);
        for G4:=1 to (t-2) do //Projede pravı sloupec
          Sx:=Sx+sqr(GPruObraz^[(I-K+G4)*width+(J+K)]-Ax);
    end;
  Sx:=sqrt(Sx/(v-1));
  result:=Sx;
end;

function Korel(r,I,J,R2,S2,width:integer;Ax,Sx:extended;ObrP2:PPoleWord):single;
var K,G1,G2,G3,G4,t:integer;
    Sy,Ay,Sumy,v:Extended;
    Korelac:single;
begin
  Sumy:=ObrP2^[(I+R2)*width+(J+S2)]; //Zjistí jas pixelu, kterı leí v okolí uzlu, tento pixel bereme z prvního obrázku v listboxu
  t:=3;v:=(2*r+1)*(2*r+1);//v je poèet pixelù v okolí uzlu sítì
  for K:=1 to r do // Cyklus pro danı poèet krunic kolem uzlu, poèet krunic= polomìr nìjvìtší krunice tj. r
    begin
      if (K>1) then t:=t+2; //Nastaví novou délku sloupce, resp. øádku
        for G1:=0 to (t-1) do //Projede spodní øádek
          Sumy:=Sumy+ObrP2^[(I+R2+K)*width+(J+S2+G1-K)];
        for G2:=0 to (t-1) do //Projede horní øádek
          Sumy:=Sumy+ObrP2^[(I+R2-K)*width+(J+S2+G2-K)];
        for G3:=1 to (t-2) do //Projede levı sloupec
          Sumy:=Sumy+ObrP2^[(I-K+G3+R2)*width+(J-K+S2)];
        for G4:=1 to (t-2) do //Projede pravı sloupec
          Sumy:=Sumy+ObrP2^[(I-K+G4+R2)*width+(J+K+S2)];
    end;
  Ay:=Sumy/v;//Aritmetickı prùmìr pro ObrP2
  t:=3;//v je poèet pixelù v okolí uzlu sítì
  Sy:=sqr(ObrP2^[(I+R2)*width+(J+S2)]-Ay);
  Korelac:=(GPruObraz^[I*width+J]-Ax)*(ObrP2^[(I+R2)*width+(J+S2)]-Ay);
  for K:=1 to r do // Cyklus pro danı poèet krunic kolem uzlu, poèet krunic= polomìr nìjvìtší krunice tj. r
    begin
      if (K>1) then t:=t+2; //Nastaví novou délku sloupce, resp. øádku
        for G1:=0 to (t-1) do //Projede spodní øádek
          begin
            Sy:=Sy+sqr(ObrP2^[(I+R2+K)*width+(J+S2+G1-K)]-Ay);
            Korelac:=Korelac+(GPruObraz^[(I+K)*width+(J-K+G1)]-Ax)*(ObrP2^[(I+R2+K)*width+(J+S2+G1-K)]-Ay);
          end;
        for G2:=0 to (t-1) do //Projede horní øádek
          begin
            Sy:=Sy+sqr(ObrP2^[(I+R2-K)*width+(J+S2+G2-K)]-Ay);
            Korelac:=Korelac+(GPruObraz^[(I-K)*width+(J-K+G2)]-Ax)*(ObrP2^[(I+R2-K)*width+(J+S2+G2-K)]-Ay);
          end;
        for G3:=1 to (t-2) do //Projede levı sloupec
          begin
            Sy:=Sy+sqr(ObrP2^[(I-K+G3+R2)*width+(J-K+S2)]-Ay);
            Korelac:=Korelac+(GPruObraz^[(I-K+G3)*width+(J-K)]-Ax)*(ObrP2^[(I-K+G3+R2)*width+(J-K+S2)]-Ay);
          end;
        for G4:=1 to (t-2) do //Projede pravı sloupec
          begin
            Sy:=Sy+sqr(ObrP2^[(I-K+G4+R2)*width+(J+K+S2)]-Ay);
            Korelac:=Korelac+(GPruObraz^[(I-K+G4)*width+(J+K)]-Ax)*(ObrP2^[(I-K+G4+R2)*width+(J+K+S2)]-Ay);
          end;
    end;
  Sy:=sqrt(Sy/(v-1));
  if (Sx=0) or (Sy=0) then Korelac:=-1
  else Korelac:=(Korelac/v)/(Sx*Sy);   //vısledná korelace mezi PruObar a ObrP2
  if(Korelac>1) then showmessage('Korelace je vìtší ne 1: '+FloatToStr(Korelac) );
  if(Korelac<-1) then showmessage('Korelace je menší ne -1: '+floattostr(Korelac)+#13#10+'....I='+inttostr(I)+', J='+inttostr(J)+', Pozice='+inttostr(I*width+J)+#13#10+'I+R2='+inttostr(I+R2)+', J+S2='+inttostr(J+S2)+', Pozice='+inttostr((I+R2)*width+J+S2) );
  result:=Korelac;
end;

//Bitmapa na prùmìr =================     1
procedure NewBitmap(width, height, n:integer; B: array of PPoleWord; Max: array of integer;Name:String);Overload;     // NewBitmap 1
var Pixels:PRGBTripleArray;
    Bitmap:TBitmap;
    Y,X,I,pru:Integer;
    prumer:Byte;
begin
  Bitmap:=CreateBmp(width, height,24);
  for Y := 0 to Height-1 do
    begin
      Pixels:=Bitmap.ScanLine[Y];
        for X := 0 to Width-1 do
          begin
            pru:=0;
            for I := 0 to n-1 do
              begin
                pru:=pru+round((B[I]^[Y*width+X]/Max[I])*255);
              end;
            prumer:=round(pru/n);
            Pixels[X].rgbtRed := prumer;
            Pixels[X].rgbtGreen := prumer;
            Pixels[X].rgbtBlue := prumer;
          end;
    end;
  Bitmap.SaveToFile(ExtractFilePath(ParamStr(0))+GNameOfFolder+'\'+Name+'.bmp');
  Bitmap.Free;
end;

procedure NewBitmap(width, height, n:integer; B: array of PPoleWord; Max: array of integer; Pr:PPoleWord;Name:string);Overload;
var Pixels:PRGBTripleArray;
    Bitmap:TBitmap;
    Y,X,I,pru:Integer;
    prumer:Byte;
begin
  Bitmap:=CreateBmp(width, height,24);
  for Y := 0 to Height-1 do
    begin
      Pixels:=Bitmap.ScanLine[Y];
        for X := 0 to Width-1 do
          begin
            pru:=0;
            for I := 0 to n-1 do
              begin
               pru:=pru+round((B[I]^[Y*width+X]/Max[I])*255);
              end;
            prumer:=round(pru/n);
            Pr^[Y*width+X]:=prumer;   // V ukazateli Pr je uloen prùmìr (ve formátu tif) z naètenıch obrazù.
            Pixels[X].rgbtRed := prumer;
            Pixels[X].rgbtGreen := prumer;
            Pixels[X].rgbtBlue := prumer;
          end;
    end;
  Bitmap.SaveToFile(ExtractFilePath(ParamStr(0))+GNameOfFolder+'\'+Name+'.bmp');
  Bitmap.Free;
end;
//Bitmapa na jeden obraz ==========     3
procedure NewBitmap(width,height:integer;max:Single;Grad:PSingleArray;Name:string);OverLoad;
var Pixels:PRGBTripleArray;
    Bitmap:TBitmap;
    J,I:integer;
    Result:Byte;
begin
  Bitmap:=CreateBmp(width, height,24);
  for J := 0 to height-1 do
          begin
            Pixels:=Bitmap.ScanLine[J];
            for I := 0 to width-1 do
              begin
                Result:=round((Grad[J*width+I]/max)*255);
                Pixels[I].rgbtBlue := Result;
                Pixels[I].rgbtGreen := Result;
                Pixels[I].rgbtRed := Result;
              end;
          end;
  Bitmap.SaveToFile(ExtractFilePath(ParamStr(0))+GNameOfFolder+'\'+Name+'.bmp');
  Bitmap.Free;
end;

//Bitmapa na jeden obraz ==========     4
procedure NewBitmap(width,height:integer;max:single;Grad:PPoleWord;Name:string);OverLoad;
var Pixels:PRGBTripleArray;
    Bitmap:TBitmap;
    J,I:integer;
    Result:Byte;
begin
  Bitmap:=CreateBmp(width, height,24);
  for J := 0 to height-1 do
          begin
            Pixels:=Bitmap.ScanLine[J];
            for I := 0 to width-1 do
              begin
                Result:=round((Grad[J*width+I]/max)*255);
                Pixels[I].rgbtBlue := Result;
                Pixels[I].rgbtGreen := Result;
                Pixels[I].rgbtRed := Result;
              end;
          end;
  Bitmap.SaveToFile(ExtractFilePath(ParamStr(0))+GNameOfFolder+'\'+Name+'.bmp');
  Bitmap.Free;
end;

procedure NewBitmapUprData(width,height:integer;max:single;Grad:PPoleWord;Name:string);OverLoad;
var Pixels:PRGBTripleArray;
    Bitmap:TBitmap;
    J,I:integer;
    Result:Byte;
begin
  Bitmap:=CreateBmp(width, height,24);
  for J := 0 to height-1 do
          begin
            Pixels:=Bitmap.ScanLine[J];
            for I := 0 to width-1 do
              begin
                Result:=round((Grad[J*width+I]/max)*255);
                Pixels[I].rgbtBlue := Result;
                Pixels[I].rgbtGreen := Result;
                Pixels[I].rgbtRed := Result;
              end;
          end;
  Bitmap.SaveToFile(ExtractFilePath(ParamStr(0))+'\'+GNameOfFolder+'\'+'UpravenáData'+'\'+Name+'.bmp');
  Bitmap.Free;
end;
               //Bitmapa pro vykreslení sítì hran
procedure NewBitmap(width,height,U:integer;max:Single;{;Brightness:PSingleArray;}{Positions:PIntegerArray}Name:string);OverLoad;  //U znaèí délku pole Positions a Brightness
var Pixels:PRGBTripleArray;
    Bitmap:TBitmap;
    J,I{K}:integer;
    K:single;
    Result:Byte;
    R:Boolean;
begin
   Bitmap:=CreateBmp(width, height,24);
   for J := 0 to height-1 do
          begin
            Pixels:=Bitmap.ScanLine[J];
            for I := 0 to width-1 do
              begin
                //R:=true;
                //for K:=0 to U-1 do
                  //begin
                    if GDicPosition.ContainsKey(J*width+I) then
                      begin
                        GDicBrightness.trygetvalue(J*width+I,K);
                        //Result:=round((Brightness[K]/max)*255);
                        Result:=round((K/max)*255);
                        Pixels[I].rgbtBlue := Result;
                        Pixels[I].rgbtGreen := Result;
                        Pixels[I].rgbtRed := Result;
                        //break;
                      end
                    else
                      begin
                       // if (R=true) then
                         // begin
                            Pixels[I].rgbtBlue := 0;
                            Pixels[I].rgbtGreen := 0;
                            Pixels[I].rgbtRed := 0;
                           // R:=false;
                         // end;
                      end;
                  //end;
              end;
          end;
  Bitmap.SaveToFile(ExtractFilePath(ParamStr(0))+GNameOfFolder+'\'+Name+'.bmp');
  Bitmap.Free;
end;

procedure TriangleNetBmp(width,height,Uk,e,p:integer; LVertex:TList);
var Pixels:PRGBTripleArray;
    Bitmap:TBitmap;
    J,I,K:integer;
    PVertex:GPPointsRecordKor;
    R:Boolean;
begin
  Bitmap:=CreateBmp(width, height,24);
  for J := 0 to height-1 do
          begin
            Pixels:=Bitmap.ScanLine[J];
            for I := 0 to width-1 do
              begin
              R:=true;
                for K:=0 to (LVertex.Count-1) do
                  begin
                    PVertex:=LVertex.Items[K];
                    if((J*width+I)=PVertex^.Position) and (K<Uk) then
                      begin
                        Pixels[I].rgbtBlue := 255;
                        Pixels[I].rgbtGreen := 255;
                        Pixels[I].rgbtRed := 255;
                        break;
                      end
                    else if((J*width+I)=PVertex^.Position) and (K>=Uk) then
                      begin
                        Pixels[I].rgbtBlue := 0;
                        Pixels[I].rgbtGreen := 255;
                        Pixels[I].rgbtRed := 0;
                        break;
                      end
                    else if (R=true) then
                      begin
                        Pixels[I].rgbtBlue := 0;
                        Pixels[I].rgbtGreen := 0;
                        Pixels[I].rgbtRed := 0;
                        R:=false;
                      end;
                  end;
              end;
          end;
  Bitmap.SaveToFile(ExtractFilePath(ParamStr(0))+GNameOfFolder+'\'+'SítTriangulace'+'_Krok_'+inttostr(GKrok)+'_Polomer_'+inttostr(GPolomer) +'_PH_'+inttostr(GPrahHodnota)+'_KrokPomSit_'+inttostr(e)+'_PolPomSit_'+inttostr(p)+'.bmp');
  Bitmap.Free;
end;

procedure NewPixel(PVertex,PPixel:GPPointsRecordKor);  //PVertex je pùvodní ukazatel, jeho hodnoty jsou zkopírovány do ukazatele PPixel
begin
  PPixel^.I:= PVertex^.I; PPixel^.J:= PVertex^.J;
  PPixel^.PosunI:= PVertex^.PosunI; PPixel^.PosunJ:=PVertex^.PosunJ;
  PPixel^.Position:=PVertex^.Position;
end;

procedure NewVertex(PVertex:GPPointsRecordKor;I,J,width:integer);
begin
  PVertex^.I:=I; PVertex^.J:=J;
  PVertex^.PosunI:=0; PVertex^.PosunJ:=0;
  PVertex^.Position:=I*width+J;
end;

procedure BitmapPosuvy(width,height:integer;Data:TList;Name:string);
var Pixels:PRGBTripleArray;
    Bitmap:TBitmap;
    J,I:integer;
    P:GPColors;
begin
  Bitmap:=CreateBmp(width, height,24);
  for J := 0 to Gheight-1 do
          begin
            Pixels:=Bitmap.ScanLine[J];
            for I := 0 to Gwidth-1 do
              begin
                P:=Data[j*Gwidth+I];
                Pixels[I].rgbtBlue := P^.B;
                Pixels[I].rgbtGreen :=P^.G;
                Pixels[I].rgbtRed := P^.R;
              end;
          end;
  Bitmap.SaveToFile(ExtractFilePath(ParamStr(0))+'\'+GNameOfFolder+'\'+'MapyPosuvù'+'\'+Name+'.bmp');
  Bitmap.Free;
end;

function Interpolation(I,J,width,w:integer;triangle:array of GPPointsRecordKor):boolean;overload;
var b,m:integer;
  Inside:boolean;
  PVertex4,PPixel:GPPointsRecordKor;
begin
  if w = 1 then
  begin
    for b := 0 to Gn-1 do
    begin
      for m := 0 to 1 do
      begin
        AGDicvertex[b].TryGetValue(Triangle[m]^.Position, Pvertex4);
        NewPixel(PVertex4,Triangle[m]);
      end;
    New(PPixel);
    Inside:=SavePixelAndCheck(I,J,width,b,Triangle, PPixel);
    end;
  end
  else if w=3 then
  begin
    for b := 0 to Gn-1 do
    begin
      for m := 0 to 2 do
      begin
        AGDicvertex[b].TryGetValue(Triangle[m]^.Position, Pvertex4);
        NewPixel(PVertex4,Triangle[m]);
      end;
    New(PPixel);
    SavePixel(I,J,width,b,Triangle, PPixel);
    end;
  end;
  result:=true;
end;

function MoveVector(arg,MaxJas,Jas:single):GPColors;
var PMoveVector:GPColors;
    x,r,g,b:single;
begin
  x:=(Jas/MaxJas)*255;
  New(PMoveVector);
  //R sloka
  if (arg<=60) or (arg>=300) then PMoveVector.R:=round(x)
  else if (arg>60) and (arg<120) then PMoveVector.R:=round(-(x/60)*arg+2*x)
  else if (arg>240) and (arg<300) then PMoveVector.R:=round((x/60)*arg-4*x)
  else PMoveVector.R:=0;
  //G sloka
  if (arg>=0) and (arg<60) then PMoveVector.G:=round((x/60)*arg)
  else if (arg>180) and (arg<240) then PMoveVector.G:=round(-(x/60)*arg+4*x)
  else if (arg>=60) and (arg<=180) then PMoveVector.G:=round(x)
  else PMoveVector.G:=0;
  //B sloka
  if (arg>120) and (arg<180) then PMoveVector.B:=round((x/60)*arg-2*x)
  else if (arg>300) and (arg<=359) then PMoveVector.B:=round(-(x/60)*arg+6*x)
  else if (arg>=180) and (arg<=300) then PMoveVector.B:=round(x)
  else PMoveVector.B:=0;
  Result:=PMoveVector;
end;
{function CheckPosition(I,J,w,width:integer;Triangle:array of GPPointsRecordKor):boolean;
var inside:boolean;
begin
  if w=1 then
  begin
    if Det2(J,I,Triangle[0]^.J,Triangle[0]^.I,Triangle[1]^.J,Triangle[1]^.I) then Inside:=Interpolation(I,J,width,w,Triangle) else Inside:=false; // ovìøíme, jestli bod leí na pøímce dané jeho dvìma nejblišími body
  end
  else if w=3 then
  begin
    Inside :=CheckPixel(I,J,Triangle); //ovìøí, jestli pixel leí v trojúhelníku
    if Inside= true then Inside:=Interpolation(I,J,width,w,Triangle); //pokud pixel neleí smau poslední vrchol, pokud pixel leí v trojúhelníku spoète jeho interpolaèní hodnotu
  end;
  result:=Inside;
end;}
function arg(re,im: single) : single;  { pocita argument komplexniho cisla }
var fi,pi : single;
begin
  pi:=3.14;
  if (abs(re)<0.000001) and (abs(im)<0.000001) then fi:=0 else
   begin
     if (re=0) then if (im>0) then fi:=90 else fi:=270;
     if (im=0) then if (re>0) then fi:= 0 else fi:=180;
     if (re>0) and (im>0) then fi:=180*arctan(im/re)/pi; //arctan vrací vısledek v radiánech, proto je nutné pøevést do stupòù
     if (re>0) and (im<0) then fi:=180*arctan(im/re)/pi+360;
     if (re<0) and (im>0) then fi:=180*arctan(im/re)/pi+180;
     if (re<0) and (im<0) then fi:=180*arctan(im/re)/pi+180;
   end;
  result:=fi;
end;

procedure SetEdit(T:Boolean;E1,E2,E3,E4,E5,E6,E7,E8:TEdit);
begin
  E1.readonly:=T; E2.readonly:=T;E3.readonly:=T;
  E4.readonly:=T; E5.readonly:=T; E6.readonly:=T;
  E7.readonly:=T;E8.readonly:=T;
end;

procedure CheckEntry(E:TEdit);
var S:integer;
begin
  if ((TrystrToInt(E.Text,S)) and (S>0)and(E.Text[1]<>'0') )=true then begin E.Color:=clGreen;E.font.Color:=clWhite; end
  else
  begin
    ShowMessage('Špatnı vstup. Vstupení parametr musí bıt pøirozené èíslo');
    E.Color:=clRed;
  end;
end;

procedure TMainForm.TvorbaSiteClick(Sender: TObject);
type PPixelRecord=^PixelRecord;
     PixelRecord = record
                    Position:Integer;
                    Brightness:Single;
                   end;
var J,I,J2,I2, FullSize,q,width, height,n,Size,Header,U,S,V:Integer;
    Grad:PSingleArray;
    f:file of byte;
    txtFile:textFile;
    Name,path:string;
    //Pro více záznamù
    A: array of Thandle;
    B: array of PPoleWord;
    B2: array of PPoleWord;
    prumer:integer;
    MaxJasPrumer, MaxJasHrany:single;     //Maximální hodnota jasu pixelu ve zprùmìrovaném obraze z obrazù hran
    //Parametry sítì
    r,e,t,H1,maxP,PrahovaHodnota,K,G1,G2,G3,G4:Integer;//r polomìr okolí bodu, e hustota sítì, maxP pozice pixelu s maximání jasem v okolí uzlu
    P:TList; //List do kterého se pøidávají záznamy o pixelech tj. záznam PixelRecord pomocí pointeru PPixelRecord
    PRecord,PRecord2:PPixelRecord;           //Ukazatel na záznam pixelu, kterı obsahuje pozici a jas pixelu
    maxB:Single;                    //Aktuální nejvyšší hodnota jasu pixelu v okolí uzlu
    PointsBrightness:PSingleArray; //Pole bodù obsahující jejich hodnotu jasu
    Ex:boolean;
begin
  TvorbaSite.Enabled:=false; Korelace.Enabled:=true;
  GDicPosition := TDictionary<integer,integer>.Create;  //vytvoøení slovníku pro pozice bodù v nich je poèítána korelace posuvu
  GDicBrightness := TDictionary<integer,single>.Create; //vytvoøení slovníku pro jasy bodù v nich je poèítána korelace posuvu
  SetEdit(true,Krok,Polomer,PrahHod,KrokPomSit,PolomerPomSit,MaxPol,PolomerKor,OkoliKorelace);
  if ((TryStrToInt(krok.Text,V)) and (TryStrToInt(Polomer.Text,V)) and (TryStrToInt(prahHod.Text,V)) and (TryStrToInt(krokPomSit.Text,V)) and (TryStrToInt(PolomerPomSit.Text,V)) and (TryStrToInt(MaxPol.Text,V)) and (TryStrToInt(PolomerKor.Text,V)) and (TryStrToInt(OkoliKorelace.Text,V)))=true  then
  begin
    if ((StrToInt(Krok.Text)>0) and (StrToInt(Polomer.Text)>0) and (StrToInt(PrahHod.Text)>0) and (StrToInt(KrokPomSit.Text)>0) and (StrToInt(PolomerPomSit.Text)>0) and (StrToInt(MaxPol.Text)>0) and (StrToInt(PolomerKor.Text)>0) and (StrToInt(OkoliKorelace.Text)>0)) =true  then
    begin
      Krok.Cursor:=CrNo; Polomer.Cursor:=CrNo;PrahHod.Cursor:=CrNo;
      KrokPomSit.Cursor:=CrNo; PolomerPomSit.Cursor:=CrNo; MaxPol.Cursor:=CrNo;
      PolomerKor.Cursor:=CrNo;OkoliKorelace.Cursor:=CrNo;
    end
    else
    begin
    SetEdit(False,Krok,Polomer,PrahHod,KrokPomSit,PolomerPomSit,MaxPol,PolomerKor,OkoliKorelace);
    showmessage('Špatné vstupní parametry. Vstupní hodnoty musí bıt pøirozená èísla!'); Exit;
    end;
  end
  else
  begin
  Krok.readonly:=False; Polomer.readonly:=False;PrahHod.readonly:=False;
  KrokPomSit.readonly:=False; PolomerPomSit.readonly:=False; MaxPol.readonly:=False;
  PolomerKor.readonly:=False;OkoliKorelace.readonly:=False;
  showmessage('Špatné vstupní parametry. Vstupní hodnoty musí bıt pøirozená èísla!'); Exit;
  end;
  path:= ExtractFilePath(ParamStr(0));
  SetCurrentDir(path);
  Ex:=true;S:=0;
  while Ex=true do
  begin
    Inc(S);
    Ex:=directoryexists('Vısledky'+inttostr(S));
  end;
  SetCurrentDir(ExtractFileDir(FileListBox1.FileName));
  GNameOfFolder:='Vısledky'+inttostr(S);
  CreateDir(ExtractFilePath(ParamStr(0))+GNameOfFolder);
  AssignFile(txtFile,path+GNameOfFolder+'\'+'ParametryInfo.txt');
  ReWrite(txtFile);
  WriteLn(txtFile,'                 Informace o nastavenıch parametrech');
  WriteLn(txtFile,'--------------------------------------------------');
  WriteLn(txtFile,'   Parametry pro tvorbu sítì');
  WriteLn(txtFile,'     Krok = '+Krok.Text);
  WriteLn(txtFile,'     Polomìr = '+polomer.Text);
  WriteLn(txtFile,'     Prahová hodnota jasu = '+PrahHod.Text);
  WriteLn(txtFile,'--------------------------------------------------');
  WriteLn(txtFile,'   Parametry pro tvorbu pomocné sítì');
  WriteLn(txtFile,'     Krok v pomocné síti = '+KrokPomSit.Text);
  WriteLn(txtFile,'     Polomìr v pomocné síti = '+PolomerPomSit.Text);
  WriteLn(txtFile,'     Maximální polomìr okolí interpolace = '+MaxPol.Text);
  WriteLn(txtFile,'--------------------------------------------------');
  WriteLn(txtFile,'   Parametry na vıpoèet korelace');
  WriteLn(txtFile,'     Polomìr v okolí uzlu = '+PolomerKor.Text);
  WriteLn(txtFile,'     Polomìr okolí na vıpoèet korelace = '+OkoliKorelace.Text);
  WriteLn(txtFile,'--------------------------------------------------');
  WriteLn(txtFile,'   Vstupní data');
  width:=0;height:=0;FullSize:=0;
  n:=0; //n znaèí poèet vybranıch obrázkù
  for J := 0 to FileListBox1.Count-1 do     //Zjistí kolik poloek z filelistboxu je vybranıch
     if FileListBox1.Selected[J] then
     begin
     if n=0 then
       begin
        if (EndsText('.tif',extractfilename(FileListBox1.FileName)))=false then
        begin
          Showmessage('Špatná vstupní data. Data musejí bıt ve formátu .tif.');
          ShellExecute(Handle, nil, PChar(Application.ExeName), nil, nil, SW_SHOWNORMAL);
          Close;
          Exit;
        end;
        AssignFile(f, FileListBox1.FileName);
        FileMode := fmOpenRead;    //Nutnı kus kódu, bez kterého by nastal I/O error 32 .... nemazat
        Reset(f);
        FullSize:=filesize(f); //2885814;
        closefile(f);
        TestTIFFile(FileListBox1.FileName,TIFFinfo);
        LoadFromFileTIF(FileListBox1.FileName,1,Images[0],TIFFinfo);
        width:=TIFFinfo[1].Width.Value;//1390;//TIFFinfo[1].Width.Value;//1390;//TIFFinfo[1].Width.Value;//1390;//TIFFinfo[1].Width.Value;//1390;//1391;
        height:=TIFFinfo[1].Height.Value;//1038;//TIFFinfo[1].Height.Value;//1038;//TIFFinfo[1].Height.Value;//1038;//1039;
       end;
     n:=n+1;
     end;
  Size:=width*height*2;
  Header:=FullSize-Size; //Pouívej pro celı obraz
  //Header:=FullSize-Size+width*838*2; //pouívej pro podmnoinu obrazu
  //height:=200; // pro podmnoinu dat, sma to pokud chceš celı obraz
  GSize:=Size;
  GHeader:=header; //pro naètení celého obrazu
  //GHeader:=header+width*838*2; // pro podmnoinu obrazu
  Gn:=n; Gwidth:=width; Gheight:=height;//height;
  // R:=Size div 2; //  R je poèet pixelù v obraze, na 1 pixel pøipadají 2 Bajty
  SetLength(A,n); SetLength(B,n); SetLength(B2,n);{SetLength(GradP,n);}SetLength(GB2,n);
  SetLength(GArrayOfNames,0);
  SetLength(GArrayOfNames,n);//Nastaví velikost polí na poèet vybranıch poloek
  q:=0;
  for I := 0 to FileListBox1.Count-1 do //Zaplní pole A handly na vybrané soubory a pole B ukazateli na pøíslušné adresy pamìti, které danım handlùm pøísluší
    begin
      if FileListBox1.Selected[I] then
        begin
          if (EndsText('.tif',extractfilename(FileListBox1.FileName)))=false then
          begin
            Showmessage('Špatná vstupní data. Data musejí bıt ve formátu .tif.');
            ShellExecute(Handle, nil, PChar(Application.ExeName), nil, nil, SW_SHOWNORMAL);
            Close;
            Exit;
          end;
          GetMem(B[q],Size); GetMem(B2[q],Size); {GetMem(GradP[q], width*height*4);} GetMem(GB2[q],Size);
          A[q]:=FileOpen(FileListBox1.Items[I],fmOpenRead);
          FileSeek(A[q],Header,0);// Z souboru pod hnadlem H pøeskoèí pøíslušnı poèet bajtù udanıch v Headeru. Èíslo 0 znaèí, e se tento poèet dat pøeskoèí od zaèátku souboru
          FileRead(A[q],B[q]^,Size);//Ze souboru pod handlem A[q] naèti poèet bajtù udanıch promìnnou Size do bufferu B2[q]^
          Move(B[q]^,B2[q]^,size);//Zkopíruje originální data do nového pole B2[q]^
          Move(B2[q]^,GB2[q]^,size);
          GArrayOfNames[q]:=extractfilename(FileListBox1.Items[I]);
          WriteLn(txtFile,    GArrayOfNames[q]);
          q:=q+1;
        end;
    end;
  CloseFile(txtFile);
    //Vypoèítání prùmìrného obrazu z obrazù hran a uloí jej do ukazatele Grad
  GetMem(Grad, width*height*4); GetMem(GPruObraz,Size);
  MaxJasPrumer:=0;                      //MaxJ = Maximální jas pixelu ve zprùmìrovaném obraze hran
  for I:=0 to (height-1) do
    begin
      for J:=0 to (width-1) do
      begin
        prumer:=0;
        for K := 0 to (n-1) do
          begin
            prumer:=prumer+B2[K]^[I*width+J];
          end;
        GPruObraz[I*width+J]:=round((prumer/n)); // prùmìrnı jas pixelu [I*width+J] ve zprùmìrovaném obraze
        if(GPruObraz[I*width+J]>MaxJasPrumer) then MaxJasPrumer:=GPruObraz[I*width+J]; //Zjistí maximální jas pixelu ve zprùmìrovaném obraze
      end;
    end;
  if (n>1) then NewBitmap(width,height,MaxJasPrumer,GPruObraz,'ZprumerovanyObraz');
    //Ve zprùmìrovaném obrazu najdu hrany
  MaxJasHrany:=Gradient(height,width,Size,GPruObraz,Grad); // Ze zprùmìrovaného obrazu v GPruObraz se napoèítají hrany a ty jsou pak uloeny do Grad
  NewBitmap(width,height,MaxJasHrany,Grad,'Hrany_ZprumerovanyObraz');   //Uloí data v Grad jako bmp
    //Tvorba sítì
  //e = krok v síti, r = polomìr okolí, t = urèuje délku øádku a sloupce v kadém okolí (nemìnit!!!), PrahovaHodnota = urèuje minimální jas pixelu, k tomu, aby byl vybrán jako prvek sítì
  e:=strtoint(Krok.Text); r:=strtoint(Polomer.Text);PrahovaHodnota:=strtoint(PrahHod.Text);
  GKrok:=e; GPolomer:=r;GPrahHodnota:=PrahovaHodnota;
  P := TList.Create;
  ProgressBar1.Position := 0;
  ProgressBar1.Min := 0;
  ProgressBar1.Max := ((height-1) div e);
  ProgressBar1.Step := round(ProgressBar1.width/((height-1) div e));
  for I:=0 to ((height-1) div e) do
    begin
    ProgressBar1.StepIt;
    ProgressBar1.Update;
    I2:=I*e;
      for J:=0 to((width-1) div e) do
        begin
          J2:=J*e;
          if (Grad^[I2*width+J2]>PrahovaHodnota) then //Pokud má uzel jas vìtší, ne je prahová hodnota
            begin
              New(PRecord);
              PRecord^.Position:=I2*width+J2;
              PRecord^.Brightness:=Grad^[I2*width+J2];
              P.Add(PRecord);
            end
          else //Nemá-li uzel jas vìtší ne prahová hodnota, hledáme v jeho okolí
            begin
              t:=3;
              for K:=1 to r do // Cyklus pro danı poèet krunic kolem uzlu, poèet krunic= polomìr nìjvìtší krunice tj. r
                begin
                    maxB:=PrahovaHodnota;
                    if (K>1) then t:=t+2; //Nastaví novou délku sloupce, resp. øádku
                    H1:=-K;
                    for G1:=0 to (t-1) do //Projede spodní øádek
                      begin
                        if(J2+H1+G1<0) or (J2+H1+G1>(width-1)) or ((I2+K)*width>(height-1)) then//Podmínka, jeslti bod v okolí uzlu leí v obrazu v platném sloupci a øádku
                        else   //Pokud bod leí ve správném okolí
                          begin
                            if(Grad^[(I2+K)*width+(J2+H1+G1)]>maxB) then
                              begin
                                maxB:= Grad^[(I2+K)*width+(J2+H1+G1)] ;  //Hodnota maximální jasu ve spodním øádku
                                maxP:= (I2+K)*width+(J2+H1+G1);          //Pozice pixelu s maximálním jasem ve spodním øádku
                              end;
                          end;
                      end;
                    for G2:=0 to (t-1) do //Projede horní øádek
                      begin
                          if(J2+H1+G2<0) or (J2+H1+G2>(width-1)) or ((I2-K)*width<0) then//Podmínka, jeslti bod v okolí uzlu leí v obrazu v platném sloupci a øádku
                        else   //Pokud bod leí ve správném okolí
                          begin
                            if(Grad^[(I2-K)*width+(J2+H1+G2)]>maxB) then
                              begin
                                maxB:= Grad^[(I2-K)*width+(J2+H1+G2)] ;  //Hodnota maximální jasu v horním øádku
                                maxP:= (I2-K)*width+(J2+H1+G2);          //Pozice pixelu s maximálním jasem v horním øádku
                              end;
                          end;
                      end;
                    for G3:=1 to (t-2) do //Projede levı sloupec
                      begin
                        if(J2-K<0) or ((I2-K+G3)*width<0) or ((I2-K+G3)*width>(height-1)) then//Podmínka, jestli bod v okolí uzlu leí v obrazu v platném sloupci a øádku
                        else      //Pokud bod leí ve správném okolí
                          begin
                            if(Grad^[(I2-K+G3)*width+(J2-K)]>maxB) then
                              begin
                                 maxB:=Grad^[(I2-K+G3)*width+(J2-K)];   //Hodnota maximální jasu v levém sloupci
                                 maxP:=(I2-K+G3)*width+(J2-K);         //Pozice pixelu s maximálním jasem v levém sloupci
                              end;
                          end;
                      end;
                    for G4:=1 to (t-2) do //Projede pravı sloupec
                      begin
                        if (J2+K>(width-1)) or ((I2-K+G4)*width<0) or ((I2-K+G4)*width>(height-1))then//Podmínka, jestli bod v okolí uzlu leí v obrazu v platném sloupci a øádku
                        else     //Pokud bod leí ve správném okolí
                          begin
                            if(Grad^[(I2-K+G4)*width+(J2+K)]>maxB)then
                              begin
                                 maxB:=Grad^[(I2-K+G4)*width+(J2+K)];   //Hodnota maximální jasu v pravém sloupci
                                 maxP:=(I2-K+G4)*width+(J2+K);         //Pozice pixelu s maximálním jasem v pravém sloupci
                              end;
                          end;
                      end;
                 if(maxB>PrahovaHodnota) then
                  begin
                    New(PRecord);
                    PRecord^.Position:=maxP;
                    PRecord^.Brightness:=maxB;
                    P.Add(PRecord);
                    break;
                  end;
                end;
            end;
        end;
    end;
    ProgressBar1.Position := 0;
  U:=P.Count;
  for I:=0 to (U-1) do
  begin
    PRecord:=P.Items[I];     //Do PRecord vloí I. záznam
    if GDicPosition.ContainsKey(PRecord^.Position)=false then
    begin
      GDicPosition.Add(PRecord^.Position,PRecord^.Position);
      GDicBrightness.Add(PRecord^.Position,PRecord^.Brightness);
    end;
  end;
    U:= GDicPosition.count; GU:=U;  //U = Poèet nalezenıch bodù, po odstranìní duplicitních záznamù
    //Uloení do bitmapy
  Name:=extractfilename(FileListBox1.FileName); //Do Name se uloí název obrázku
  delete(Name,Length(Name)-3,4);  //Z Name se oddstraní poslední 4 znaky tj. .tif
  if (n=1) then  NewBitmap(width,height,U,MaxJasPrumer{PointsBrightness}{GPointsPosition},'SítZHran_Krok_'+inttostr(e)+'_polomer_'+inttostr(r)+'_PH_'+inttostr(PrahovaHodnota)+'_'+Name) //Uloí sí do bitmapy a pøipojí název obrázku
  else NewBitmap(width,height,U,MaxJasHrany{,PointsBrightness]{GPointsPosition},'SítZHran_Krok_'+inttostr(e)+'_polomer_'+inttostr(r)+'_PH_'+inttostr(PrahovaHodnota)); //Uloí sí ze zprùmìrovaného obrázku do bitmapy
    //Uvolnìní pamìti
  FreeMem(Grad); GDicBrightness.Clear; GDicBrightness.Destroy; //FreeMem(PointsBrightness,U*4); //FreeMem(PointsPosition); //Uvolní zabranou pamì
  for I := 0 to (U - 1) do
  begin
    PRecord := P.Items[I];
    Dispose(PRecord);   //Uvolní jednotlivé záznamy
  end;
  P.Free;   //Uvolní list P
  for I:=0 to (n-1) do//Uvolní jednotlivé buffery, n je poèet vybranıch obrázkù
  begin
    FreeMem(B[I]);FreeMem(B2[I]);
    fileclose(A[I]);
  end;
  setlength(A,0); setlength(B,0); setlength(B2,0);
  showmessage('Sí dokonèena.');
 // ReportMemoryLeaksOnShutdown := True;
end;

procedure TMainForm.KorelaceClick(Sender: TObject);
type PPointsRecord=^PointsRecord;
     PointsRecord = record I,J:Integer; end;
var J,I,J2,I2,L,FullSize,q,width,height,U,Header,Size, BodI, BodJ,b,G:Integer;
    Name:string;
    Up:boolean;
    LUzly,Lbody:Tlist;
    PRecordU,PRecordB:PPointsRecord;
    MaxKorelace,KorelaceR:single;
    ObrP2 : PPoleWord; //Obraz = Data prvního obrazu v listboxu, ke kterému poèítáme korelaci, ObrP = Obraz
    //Parametry sítì
    Ax,Sx:extended;//Ax=aritmetickı prùmìr hodnot pixelù leících v okolí uzlu, Sx=smìrodatná odchylka hodnot pixelù v okolí uzlu sítì
    r,Rk,a,e,t,H1,K,G1,G2,G3,G4:Integer;//r polomìr okolí bodu, e hustota sítì,
    PRecordKor,PRecordKor2:GPPointsRecordKor;
 begin
  Korelace.Enabled:=false; Interpolace.Enabled:=true;
  width:=Gwidth; height:=Gheight; Header:=Gheader; U:=GU; Size:=GSize;
  SetLength(GPointsKor, Gn);//nastaví délku pole pro TListy na poèet vybranıch obrázkù Gn
  SetLength(AGDicVertex, Gn);//nastaví délku pole pro slovníky na poèet vybranıch obrázkù Gn
  ProgressBar1.Position := 0;
  ProgressBar1.Update;
  ProgressBar1.Min := 0;
  ProgressBar1.Max :=Gn ;
  ProgressBar1.Step := 1;//round(ProgressBar1.width/height);
  for b := 0 to Gn-1 do
  begin
  //ProgressBar1.Position := 0;
  ProgressBar1.StepIt;
  ProgressBar1.Update;
  LUzly:=TList.Create; LBody:=Tlist.Create;
  GPointsKor[b]:=TList.Create;
  AGDicVertex[b] := TDictionary<integer, GPPointsRecordKor>.Create;  //vytvoøení slovníku pro vrcholy interpolaèní sítì
  GetMem(ObrP2,Size);
  Move(GB2[b]^,ObrP2^,size);//Zkopíruje originální data do nového bufferu ObrP2^, kterı v sobì obsahuje data obrázku q-tého
     //Projede body v obrazu ObrP2 a v bodech, jejich souøadnice odpovídají bodùm v síti spoète pro kadı bod jejich okolí korelaci
  r:=strtoint(PolomerKor.Text);a:=strtoint(OkoliKorelace.Text); Rk:=a; // r = Polomìr okolí uzlu sítì, Rk = polomìr okolí uzlu sítì, v nìm hledáme body na vıpoèet korelace, Up pokud leí pixel s maximální korelací na hranici okolí v nìm hledáme korelaci, pak toto okolí zvyš o +1
  G:=0;
  for I:=0 to (height-1) do
    begin
      //ProgressBar1.StepIt;
      //ProgressBar1.Update;
      for J:=0 to (width-1) do
      begin
         if ( GDicPosition.ContainsKey(I*width+J) {I*width+J=GPointsPosition[L]}) And ( ( ((I-r)>=0) And ((I+r)<=(height-1)) ) And ( ((J-r)>=0) And ((J+r)<=(width-1)) ) ) then   //Zjistí, jestli je bod uzlem sítì a jestli jeho celé okolí leí v obraze
          begin
            Ax:=AritmetickyPrumer(r,I,J,width);
            Sx:=SmerodatnaOdchylka(r,I,J,width,Ax);
            K:=0; t:=1;
            while ((K<=Rk) OR (Up=True) ) do
              begin
                if K=0 then
                  begin
                    KorelaceR:=Korel(r,I,J,K,K,width,Ax,Sx,ObrP2);
                    MaxKorelace:=KorelaceR;
                    BodI:=I; BodJ:=J;
                  end
                else
                begin
                t:=t+2; //Nastaví novou délku sloupce, resp. øádku
                H1:=-K; Up:=false;  //Up pokud leí pixel s maximální korelací na hranici okolí v nìm hledáme korelaci, pak toto okolí zvyš o +1
                for G1:=0 to (t-1) do //Projede spodní øádek
                  begin
                    if (J+H1+G1-r<0) or (J+H1+G1+r>(width-1)) or ((I+K+r)>(height-1)) or ((I+K-r)<0) then // Podmínka, jestli okolí bodu v nìm poèítáme korelaci v okolí uzlu leí v obraze. poslední ètyøi podmínky ovìøují jestli lze spoèíst korelaci s polomìrem r k uzlu
                    else
                      begin
                        KorelaceR:=Korel(r,I,J,K,(-K+G1),width,Ax,Sx,ObrP2);
                        if (korelaceR>MaxKorelace) then
                          begin
                            MaxKorelace:=KorelaceR;
                            BodI:=I+K; BodJ:=J+H1+G1;
                            if K>=Rk then Up:=True; // pokud bod s maximální korelací, leí v na hranici okolí uzlu, pak toto oklí zvıšíme o +1 (viz. øádek níe)
                          end;
                      end;
                    end;
                for G2:=0 to (t-1) do //Projede horní øádek
                  begin
                    if (J+H1+G2-r<0) or (J+H1+G2+r>(width-1)) or ((I-K+r)>(height-1)) or ((I-K-r)<0) then // Podmínka, jestli okolí bodu v nìm poèítáme korelaci v okolí uzlu leí v obraze.  poslední ètyøi podmínky ovìøují jestli lze spoèíst korelaci s polomìrem r k uzlu
                    else
                      begin
                        KorelaceR:=Korel(r,I,J,-K,(-K+G2),width,Ax,Sx,ObrP2);
                        if (korelaceR>MaxKorelace) then
                          begin
                            MaxKorelace:=KorelaceR;
                            BodI:=I-K;BodJ:=J+H1+G2;
                            if K>=Rk then Up:=True; // pokud bod s maximální korelací, leí v na hranici okolí uzlu, pak toto oklí zvıšíme o +1 (viz. øádek níe)
                          end;
                      end;
                    end;
                for G3:=1 to (t-2) do //Projede levı sloupec
                    begin
                    if (J-K-r<0) or ((I-K+G3+r)>(height-1)) or ((I-K+G3-r)<0) or (J-K+r>(width-1)) then // Podmínka, jestli okolí bodu v nìm poèítáme korelaci v okolí uzlu leí v obraze.  poslední ètyøi podmínky ovìøují jestli lze spoèíst korelaci s polomìrem r k uzlu
                    else
                      begin
                        KorelaceR:=Korel(r,I,J,(-K+G3),-K,width,Ax,Sx,ObrP2);
                        if (korelaceR>MaxKorelace) then
                          begin
                            MaxKorelace:=KorelaceR;
                            BodI:=I-K+G3; BodJ:=J-K;
                            if K>=Rk then Up:=True; // pokud bod s maximální korelací, leí v na hranici okolí uzlu, pak toto oklí zvıšíme o +1 (viz. øádek níe)
                          end;
                      end;
                    end;
                for G4:=1 to (t-2) do //Projede pravı sloupec
                    begin
                    if (J+K+r>(width-1)) or (J+K-r<0) or ((I-K+G4-r)<0) or ((I-K+G4+r)>(height-1)) then  // Podmínka, jestli okolí bodu v nìm poèítáme korelaci v okolí uzlu leí v obraze.  poslední ètyøi podmínky ovìøují jestli lze spoèíst korelaci s polomìrem r k uzlu
                    else
                      begin
                        KorelaceR:=Korel(r,I,J,(-K+G4),K,width,Ax,Sx,ObrP2);
                        if (korelaceR>MaxKorelace) then
                          begin
                            MaxKorelace:=KorelaceR;
                            BodI:=I-K+G4; BodJ:=J+K;
                            if K>=Rk then Up:=True; // pokud bod s maximální korelací, leí v na hranici okolí uzlu, pak toto oklí zvıšíme o +1 (viz. øádek níe)
                          end;
                      end;
                    end;
                end;
                inc(K);
              end;
          //Naètení pozice uzlu v nìm je poèítána korelace s bodem leícím v okolí tohoto uzlu do listù záznamù
          Rk:=a;
          New(PRecordU); New(PRecordB);
          PRecordU^.I:=I; PRecordU^.J:=J;
          PRecordB^.I:=BodI; PRecordB^.J:=BodJ;
          LUzly.Add(PRecordU); LBody.Add(PRecordB);
          end;
      end;
    end;
  U:=LUzly.Count; //U= poèet bodù v nich byla spoètena korelace
  //NAÈTENÍ DAT DO SLOVNÍKU AGDicCVertex a listu GPointsKor
  for I:=0 to (U-1) do
  begin
      PRecordU:=LUzly.Items[I]; PRecordB:=LBody.Items[I];
      New(PRecordKor);
      PRecordKor^.I:=PRecordU^.I; PRecordKor^.J:=PRecordU^.J;
      PRecordKor^.PosunI:=-PRecordU^.I+PRecordB^.I;
      PRecordKor^.PosunJ:=-PRecordU^.J+PRecordB^.J;
      PrecordKor^.Position:=PRecordU^.I*width+PRecordU^.J;
      GPointsKor[b].Add(PRecordKor);
      New(PRecordKor2);
      PRecordKor2^.I:=PRecordU^.I; PRecordKor2^.J:=PRecordU^.J;
      PRecordKor2^.PosunI:=-PRecordU^.I+PRecordB^.I;
      PRecordKor2^.PosunJ:=-PRecordU^.J+PRecordB^.J;
      PrecordKor2^.Position:=PRecordU^.I*width+PRecordU^.J;
      if AGDicVertex[b].ContainsKey(PrecordKor2^.Position)=false then
        AGDicVertex[b].Add(PrecordKor2^.Position, PrecordKor2); // nahrání pùvodních vrcholù do slovníku interpolaèní sítì
  end;
  for I := 0 to U-1 do  //Uvolnìní pamìti vyhrazené na pointery v listu
  begin
    Dispose(LUzly.Items[I]);
    Dispose(LBody.Items[I]);
  end;
  LUzly.Free; LBody.Free;//smae listy
  //Uvolnìní pamìti
  FreeMem(ObrP2); //Uvolní zabranou pamì, GPointsPosition - list uzlù sítì
  end;
  showmessage('Korelace dokonèena.');
  GDicPosition.Clear; GDicPosition.destroy;
  //FreeMem(GPointsPosition,GU*4);
 // ReportMemoryLeaksOnShutdown := True;
end;

procedure TMainForm.InterpolaceClick(Sender: TObject);
var width, height,I,J,I2,J2,m,e,r,t,t2,L,H1,K,K2,Uk,G1,G2,G3,G4,w,P,Q,C,b,Key,S:Integer; //e = krok v pomocné síti, r = polomìr okolí bodu v pomocné síti
    LVertex,LFakeVertex:Tlist;//LVertex = List vrcholù trojúhelníkové sítì,
    FakeVertex:GPFVertex;//ukazatel na vrchol, kterı se nehodí k aktuální triangulaci
    PPixel, PVertex,PVertex2,A,A2,A3:GPPointsRecordKor; // A,B,C pøedtsavují vrcholy trojúhelníku, kterı pouíváme pro interpolaci hodnot pixelu v trojúhelníku ABC
    Clear,Inside,FirstK,Up:Boolean;
    Triangle:array[0..2] of GPPointsRecordKor;
    D:array of array of GPPointsRecordKor;
    MaxPolomer:byte;
Label StartNewK;
begin
  Interpolace.Enabled:=false; UlozPosuvy.Enabled:=true; Posuvy.Enabled:=true;  BilinearniInterpolace.Enabled:=true;
  width:=Gwidth; height:=Gheight; LVertex:=Tlist.Create;
  Uk:= GPointsKor[0].Count; // Uk = Poèet záznamù v listu bodù, ve kterıch poèítáme korelaci
  SetLength(GLPixelData,Gn);
  MaxPolomer:=StrToInt(MaxPol.Text);
  for I:=0 to (Uk-1) do //zkopíruje hodnoty v listu GPointsKor do nového lsitu vrhcolù trojúhelníkù trojúhelníkové sítì LVertex
  begin
    New(PVertex2);
    NewPixel(GPointsKor[0].Items[I],PVertex2);
    LVertex.Add(PVertex2);
  end;
  ProgressBar1.Position := 0;
  ProgressBar1.Update;
  ProgressBar1.Min := 0;
  //Do LVertex pøidáme pomocnou sí bodù
  e:=strtoint(KrokPomSit.Text);r:=strtoint(PolomerPomSit.Text);   //e = krok v pomocné síti, r = polomìr okolí bodu v pomocné síti
  ProgressBar1.Step := round(ProgressBar1.width/(((height-1) div e)+height));
  ProgressBar1.Max :=((height-1) div e)+ height;
  for I:=0 to ((height-1) div e) do
  begin
  ProgressBar1.StepIt;
  ProgressBar1.Update;
  I2:=I*e;
    for J:=0 to ((width-1) div e) do
    begin
      J2:=J*e;
      if (I=0) or (J=0) then  //pokud pixel leí v prvním øádku, nebo prvním sloupci je automaticky uloen
      begin
        New(PVertex2);
        NewVertex(PVertex2,I2,J2,width); //Procedura, která do ukazatele pvertex nahraje pozici bodu pomocného uzlu, jeho posunutí a souøadnice I,J
        for b := 0 to Gn-1 do
        begin
          New(PVertex); NewVertex(PVertex,I2,J2,width);
          if b=0 then LVertex.Add(PVertex2); //pro sít tvoøenou uzly sítì
          if AGDicVertex[b].ContainsKey(PVertex^.Position)=false then AGDicVertex[b].Add(PVertex^.Position, PVertex); //do slovníku vrcholù interpolaèní sítì pøidáme body pomocné sítì
        end;
       continue;
      end;
      Clear:=HelpVertex(Uk,width,I2,J2,0,0,0);
      if clear=false then continue// pokud je uzel pomocné sítì identickı s uzlem pùvodní sítì, pak se pøejde na další uzel
      else //Není-li pomocnı uzel identickı s ádnım uzlem pomocné sítì, pak prohledáme okolí pomocného uzlu, abychom ovìøili, e v nìm není ádnı bod, kterı by byl uzlem pùvodní sítì
            begin
              t:=3;
              for K:=1 to r do // Cyklus pro danı poèet krunic kolem uzlu, poèet krunic= polomìr nìjvìtší krunice tj. r
                begin
                    if (K>1) then t:=t+2; //Nastaví novou délku sloupce, resp. øádku
                    H1:=-K;
                    for G1:=0 to (t-1) do //Projede spodní øádek
                    begin
                      if (I2+K>=0) and (I2+K<=height-1) and (J2-K+G1>=0) and (J2-K+G1<=width-1) then
                        begin
                          Clear:=HelpVertex(Uk,width,I2,J2,K,-K+G1,0);
                          if clear=false then break;
                        end;
                    end;
                    if clear=false then break;
                    for G2:=0 to (t-1) do //Projede horní øádek
                    begin
                      if (I2-K>=0) and (I2-K<=height-1) and (J2-K+G2>=0) and (J2-K+G2<=width-1) then
                        begin
                          Clear:=HelpVertex(Uk,width,I2,J2,-K,-K+G2,0);
                          if clear=false then break;
                        end;
                    end;
                    if clear=false then break;
                    for G3:=1 to (t-2) do //Projede levı sloupec
                    begin
                      if (I2-K+G3>=0) and (I2-K+G3<=height-1) and (J2-K>=0) and (J2-K<=width-1) then
                      begin
                        Clear:=HelpVertex(Uk,width,I2,J2,-K+G3,-K,0);
                        if clear=false then break;
                      end;
                    end;
                    if clear=false then break;
                    for G4:=1 to (t-2) do //Projede pravı sloupec
                    begin
                      if (I2-K+G4>=0) and (I2-K+G4<=height-1) and (J2+K>=0) and (J2+K<=width-1) then
                      begin
                        Clear:=HelpVertex(Uk,width,I2,J2,-K+G4,K,0);
                        if clear=false then break;
                      end;
                    end;
                    if clear=false then break;
                end;
                if clear=true then //pokud v okolí uzlu pomocné sítì neleí ádnı bod pùvodní sítì, pak je tento bod pøidán do listu vrcholù trojùhelníkové sítì LVertex
                begin
                  New(PVertex2);
                  NewVertex(PVertex2,I2,J2,width);//Procedura, která do ukazatele pvertex nahraje pozici bodu pomocného uzlu, jeho posunutí a souøadnice I,J
                  for b := 0 to Gn-1 do
                  begin
                    New(PVertex); NewVertex(PVertex,I2,J2,width);
                    if b=0 then LVertex.Add(PVertex2);//pro sít tvoøenou uzly sítì
                    if AGDicVertex[b].ContainsKey(PVertex^.Position)=false then AGDicVertex[b].Add(PVertex^.Position, PVertex); //do slovníku vrcholù interpolaèní sítì pøidáme body pomocné sítì
                  end;
                end;
            end;
    end;
  end;
  //pøidáme do pomocné sítì body na hranici obrazu
  //poslední sloupec
  for I:=0 to ((height-1) div e) do
  begin
    I2:=I*e;
    for b := 0 to Gn-1 do
    begin
      if AGDicVertex[b].ContainsKey((I2*width)+(width-1))=false then
      begin
        if b=0 then
        begin
          New(PVertex2);
          NewVertex(PVertex2,I2,width-1,width);
          if b=0 then LVertex.Add(PVertex2);  //pro sít tvoøenou uzly sítì
        end;
      New(PVertex);
      NewVertex(PVertex,I2,width-1,width);
      AGDicVertex[b].Add(PVertex^.Position, PVertex); //do slovníku vrcholù interpolaèní sítì pøidáme body pomocné sítì
      end;
    end;
  end;
  //poslední øádek
  for I:=0 to ((width-1) div e) do
  begin
    I2:=I*e;
    for b := 0 to Gn-1 do
    begin
      if AGDicVertex[b].ContainsKey((height-1)*width+I2)=false then
      begin
        if b=0 then
        begin
          New(PVertex2);
          NewVertex(PVertex2,height-1,I2,width);
          LVertex.Add(PVertex2); //pro sít tvoøenou uzly sítì
        end;
        New(PVertex);
        NewVertex(PVertex,height-1,I2,width);
        AGDicVertex[b].Add(PVertex^.Position, PVertex); //do slovníku vrcholù interpolaèní sítì pøidáme body pomocné sítì
      end;
    end;
  end;
  //pravı dolní roh
  for b := 0 to Gn-1 do
  begin
    if AGDicVertex[b].ContainsKey((height-1)*width+(width-1))=false then
    begin
      if b=0 then
      begin
        New(PVertex2);
        NewVertex(PVertex2,height-1,width-1,width);
        LVertex.Add(PVertex2);     //pro sít tvoøenou uzly sítì
      end;
    New(PVertex); NewVertex(PVertex,height-1,width-1,width);
    AGDicVertex[b].Add(PVertex^.Position, PVertex); //do slovníku vrcholù interpolaèní sítì pøidáme body pomocné sítì
    end;
  end;
  TriangleNetBmp(width,height,Uk,e,r,LVertex);//LVertex);    //Uloení triangulaèní sítì
  //UVOLNÌNÍ PAMÌTI
  for I:=0 to (LVertex.Count-1) do
  begin
    Dispose(LVertex.Items[I]);             //smae jednotlivé pointery na vrcholy trojúhelníkové sítì
  end;
  Lvertex.Free;  //smae list vrcholù trojúhelníkové sítì
    //=======================================Samotná interpolace dat pomocí interpolaèních rovin=============================================

 // ProgressBar1.Position := 0;
 // ProgressBar1.Update;
  //ProgressBar1.Min := 0;
 // ProgressBar1.Max :=height ;
  ProgressBar1.Step := 1;//round(ProgressBar1.width/height);
  New(A);New(A2);New(A3);
  for b := 0 to Gn-1 do
  begin
    GLPixelData[b]:=Tlist.Create;//vytvoøí list na pixely v obraze
  end;
  for I:=0 to (height-1) do
  begin
    ProgressBar1.StepIt;
    ProgressBar1.Update;
    for J:=0 to (width-1) do
    begin
      LFakeVertex:=TList.Create; //vytvoøí list pro nevhodné vrcholy
      if AGDicvertex[0].ContainsKey(I*width+J)=true then
      begin
        for b := 0 to Gn-1 do
        begin
          AGDicvertex[b].TryGetValue(I*width+J,PVertex);
          New(PPixel);
          NewPixel(PVertex,PPixel);//hodnoty z pointeru PVertex jsou zkopírovány do pointeru PPixel
          GLPixeldata[b].Add(PPixel);//Pokud je nìkterı pixel totonı s vrcholem interpolaèní sítì, pak se tento vrchol nahraje do listu LPixelData
        end;
      end
      else
      begin
        FirstK:=true;
        StartNewK:
        if FirstK = true then
        begin
          K2:=1; t:=3; w:=0;
          A^.I:=-1;  A^.J:=-1; A^.Position:=-1; A^.PosunI:=-1; A^.PosunJ:=-1;
          A2^.I:=-1;  A2^.J:=-1; A2^.Position:=-1; A2^.PosunI:=-1; A2^.PosunJ:=-1;
          A3^.I:=-1;  A3^.J:=-1; A3^.Position:=-1; A3^.PosunI:=-1; A3^.PosunJ:=-1;
          Triangle[0]:=A; Triangle[1]:=A2; Triangle[2]:=A3;
        end
        else
        begin
          w:=1; t:=t2;
        end;
        //===================================================================PROHLEDÁVÁNÍ OKOLÍ=====================================
        for K:=K2 to 500 do // Cyklus pro danı poèet krunic kolem uzlu, poèet krunic= polomìr nìjvìtší krunice tj. r
          begin
            if (K>1) then t:=t+2; //Nastaví novou délku sloupce, resp. øádku
            if K>MaxPolomer{4.5*e} then //pokud se prohledalo pøíliš velké okolí, pak sma nalezené 2 nejbliší body, 2. nejbliší bod ulo jako FakeVertex (tj. vrchol nevhodnı na triangulaci)
            begin
              New(FakeVertex);
              FakeVertex^.I:=Triangle[1]^.I; FakeVertex^.J:=Triangle[1]^.J;
              LFakeVertex.add(FakeVertex);
              FirstK:=false;
              goto StartNewK;
            end;
            Inside:=false;
            for G1:=0 to (t-1) do //Projede spodní øádek
            begin
              if (I+K>=0) and (I+K<=height-1) and (J-K+G1<=width-1) and (J-K+G1>=0) then//zjistí, jestli danı pixel leí v obraze
              begin
              if (AGDicvertex[0].ContainsKey((I+K)*width+J-K+G1)=true) and (CheckFakeVertex(LFakeVertex,I,J,K,-K+G1)=false) then
              begin
                if (Triangle[0]^.I<>I+K) or (Triangle[0]^.J<>J-K+G1)  then
                begin
                  AGDicvertex[0].TryGetValue((I+K)*width+J-K+G1,PVertex);
                  NewPixel(PVertex,Triangle[w]);
                  if w=1 then
                  begin
                   if Det2(J,I,Triangle[0]^.J,Triangle[0]^.I,Triangle[1]^.J,Triangle[1]^.I) then // ovìøíme, jestli bod leí na pøímce dané jeho dvìma nejblišími body
                   begin
                    Inside:=Interpolation(I,J,width,w,Triangle);
                    break;
                   end;
                   K2:=K;  t2:=t;
                  end;
                  inc(w);
                  if w = 3 then  //pokud w =3, pak jsme našli 3 nejbliší body a musíme dále ovìøit, jestli náš pixel, leí uvnitø trojúhelníku tvoøeného tìmito body
                  begin
                    Inside :=CheckPixel(I,J,Triangle); //ovìøí, jestli pixel leí v trojúhelníku
                    if Inside= false then Dec(w)  //pokud neleí smau poslední vrchol
                    else                    //pokud pixel leí v trojúhelníku spoète jeho interpolaèní hodnotu
                    begin
                      Inside:=Interpolation(I,J,width,w,Triangle);
                      break;
                    end;
                  end;
                end;//konec bloku na zjišování, zda bod leí v listu vrcholù triangulaèní sítì
              end;
              end;//konec bloku, zda pixel leí v obraze
            end; //konec procházení øádku
            if Inside = true then break;
            for G2:=0 to (t-1) do //Projede horní øádek
            begin
              if (I-K>=0) and (I-K<=height-1) and (J-K+G2<=width-1) and (J-K+G2>=0) then//zjistí, jestli danı pixel leí v obraze
              begin
              if (AGDicvertex[0].ContainsKey((I-K)*width+J-K+G2)=true) and (CheckFakeVertex(LFakeVertex,I,J,-K,-K+G2)=false) then
              begin
                if (Triangle[0]^.I<>I-K) or (Triangle[0]^.J<>J-K+G2) then
                begin
                  AGDicvertex[0].TryGetValue((I-K)*width+J-K+G2,PVertex);
                  NewPixel(PVertex,Triangle[w]);
                  if w = 1 then // pokud jsme našli dva nejbliší uzly, ovìøíme, jestli náš bod, pro kterı poèítáme interpolaci leí na úseèce mezi tìmito body
                  begin
                   if Det2(J,I,Triangle[0]^.J,Triangle[0]^.I,Triangle[1]^.J,Triangle[1]^.I) then // ovìøíme, jestli bod leí na pøímce dané jeho dvìma nejblišími body
                   begin
                      Inside:=Interpolation(I,J,width,w,Triangle);
                      break;
                   end;
                   K2:=K; t2:=t;
                  end;
                  inc(w);
                  if w = 3 then  //pokud w =3, pak jsme našli 3 nejbliší body a musíme dále ovìøit, jestli náš pixel, leí uvnitø trojúhelníku tvoøeného tìmito body
                  begin
                    Inside :=CheckPixel(I,J,Triangle); //ovìøí, jestli pixel leí v trojúhelníku
                    if Inside= false then Dec(w)  //pokud neleí smau poslední vrchol
                    else                    //pokud pixel leí v trojúhelníku spoète jeho interpolaèní hodnotu
                    begin
                      Inside:=Interpolation(I,J,width,w,Triangle);
                      break;
                    end;
                  end;
                end;
              end;
              end;
            end;
            if Inside = true then break;
            for G3:=1 to (t-2) do //Projede levı sloupec
            begin
              if (I-K+G3>=0) and (I-K+G3<=height-1) and (J-K<=width-1) and (J-K>=0) then//zjistí, jestli danı pixel leí v obraze
              begin
                if (AGDicvertex[0].ContainsKey((I-K+G3)*width+J-K)=true) and (CheckFakeVertex(LFakeVertex,I,J,-K+G3,-K)=false) then
                begin
                if (Triangle[0]^.I<>I-K+G3) or (Triangle[0]^.J<>J-K) then
                begin
                  AGDicvertex[0].TryGetValue((I-K+G3)*width+J-K,PVertex);
                  NewPixel(PVertex,Triangle[w]);
                  if w = 1 then // pokud jsme našli dva nejbliší uzly, ovìøíme, jestli náš bod, pro kterı poèítáme interpolaci leí na úseèce mezi tìmito body
                  begin
                   if Det2(J,I,Triangle[0]^.J,Triangle[0]^.I,Triangle[1]^.J,Triangle[1]^.I) then // ovìøíme, jestli bod leí na pøímce dané jeho dvìma nejblišími body
                   begin
                    Inside:=Interpolation(I,J,width,w,Triangle);
                    break;
                   end;
                   K2:=K; t2:=t;
                  end;
                  inc(w);
                  if w = 3 then  //pokud w =3, pak jsme našli 3 nejbliší body a musíme dále ovìøit, jestli náš pixel, leí uvnitø trojúhelníku tvoøeného tìmito body
                  begin
                    Inside :=CheckPixel(I,J,Triangle); //ovìøí, jestli pixel leí v trojúhelníku
                    if Inside= false then Dec(w)  //pokud neleí smau poslední vrchol
                    else                    //pokud pixel leí v trojúhelníku spoète jeho interpolaèní hodnotu
                    begin
                      Inside:=Interpolation(I,J,width,w,Triangle);
                      break;
                    end;
                  end;
                end;
                end;
              end;
            end;
            if Inside = true then break;
            for G4:=1 to (t-2) do //Projede pravı sloupec
            begin
              if (I-K+G4>=0) and (I-K+G4<=height-1) and (J+K<=width-1) and (J+K>=0) then//zjistí, jestli danı pixel leí v obraze
              begin
                if (AGDicvertex[0].ContainsKey((I-K+G4)*width+J+K)=true) and (CheckFakeVertex(LFakeVertex,I,J,-K+G4,K)=false) then
                begin
                if (Triangle[0]^.I<>I-K+G4) or (Triangle[0]^.J<>J+K) then
                begin
                  AGDicvertex[0].TryGetValue((I-K+G4)*width+J+K,PVertex);
                  NewPixel(PVertex,Triangle[w]);
                  if w = 1 then // pokud jsme našli dva nejbliší uzly, ovìøíme, jestli náš bod, pro kterı poèítáme interpolaci leí na úseèce mezi tìmito body
                  begin
                   if Det2(J,I,Triangle[0]^.J,Triangle[0]^.I,Triangle[1]^.J,Triangle[1]^.I) then // ovìøíme, jestli bod leí na pøímce dané jeho dvìma nejblišími body
                   begin
                    Inside:=Interpolation(I,J,width,w,Triangle);
                    break;
                   end;
                   K2:=K; t2:=t;
                  end;
                  inc(w);
                  if w = 3 then  //pokud w =3, pak jsme našli 3 nejbliší body a musíme dále ovìøit, jestli náš pixel, leí uvnitø trojúhelníku tvoøeného tìmito body
                  begin
                    Inside :=CheckPixel(I,J,Triangle); //ovìøí, jestli pixel leí v trojúhelníku
                    if Inside= false then Dec(w)  //pokud neleí smau poslední vrchol
                    else                    //pokud pixel leí v trojúhelníku spoète jeho interpolaèní hodnotu
                    begin
                      Inside:=Interpolation(I,J,width,w,Triangle);
                      break;
                    end;
                  end;
                end;
                end;
              end;
            end;
            if Inside = true then break;
          end;//konec K cyklu, cyklus na prohledávání okolí
      end;
    for c:=0 to (LFakeVertex.Count-1) do Dispose(LFakeVertex.Items[C]);  //smae pointery na nevhodné vrcholy pro aktuální pixel
    LFakeVertex.Free;      //smae list nevhodnıch pixelù
    end; //Konec J cyklu
  end; //Konec I cyklu
  Dispose(A); Dispose(A2);Dispose(A3);
  SetLength(D, Gn); S:=0;
  for I := 0 to Gn-1 do setLength(D[I],AGDicVertex[I].Count);
  Up:=false;
  for I := 0 to Height-1 do
  begin
    for J := 0 to width-1 do
    begin
      for b := 0 to Gn-1 do
      begin
        if AGDicVertex[b].ContainsKey(I*width+J)=true then
          begin
            Up:=true;
            AGDicVertex[b].TryGetValue(I*width+J,Pvertex);
            D[b,S]:=Pvertex;
          end;
      end;
      if Up=true then inc(S);
      Up:=false;
    end;
  end;
  for b := 0 to Gn-1 do           //UVOLNÌNÍ PAMÌTI ALOKOVANÉ NA SLOVNÍKY VRCHOLÙ INTERPOLAÈNÍ SÍTÌ
  begin
    AGDicvertex[b].Clear; //AGDicvertex[b].Free;
    AGDicvertex[b].Destroy;
  end;
  for b := 0 to Gn-1 do
  begin
    for S := 0 to Length(D[b])-1 do dispose(D[b,S]);
  end;
  for b := 0 to Gn-1 do setlength(D[b],0);
  setlength(D,0);

  for I := 0 to Gn-1 do
  begin
    for J := 0 to (Uk-1) do dispose(GPointsKor[I].Items[J]);
    GPointsKor[I].Free;
  end;
  ProgressBar1.Position := 0;
  SetLength(GPointsKor,0);
  Setlength(AGDicVertex,0);
  showmessage('Interpolace dat dokonèena.');
 // ReportMemoryLeaksOnShutdown := True;
end;

procedure TMainForm.PosuvyClick(Sender: TObject);
var I,b:integer;
    PPixel:GPPointsRecordKor;
    Jas,ar:Single;
    MaxJas:Array of single;
    Data:Tlist;
    statFile:TextFile;
    str:string;
    ArPruJ,ArPruI,ArPruN,RozptylJ,RozptylI,SmerOdI,SmerOdJ:array of extended;
begin
  Posuvy.Enabled:=false; MeritkoMapyPosuvu.enabled:=true;
  //Vypíše statistické vyhodnocení posuvù pro jednotlivé obrázky
  AssignFile(statFile,ExtractFilePath(ParamStr(0))+GNameOfFolder+'\'+'StatistickéVyhodnocení.txt');
  ReWrite(statFile);
  CreateDir(ExtractFilePath(ParamStr(0))+'\'+GNameOfFolder+'\'+'MapyPosuvù');
  WriteLn(statFile,';Ar. prùmìr posuvù ve vodorovném smìru;Ar. prùmìr posuvù ve svislém smìru;'+'Ar. prùmìr normy posuvù;Rozptyl posuvù ve vodorovném smìru;Rozptyl posuvù ve svislém smìru;Smìrodatná odchylka posuvù ve vodorovném smìru;Smìrodatná odchylka posuvù ve svislém smìru;Nejvìtší posuv ve smyslu normy');
  SetLength(MaxJas,Gn); SetLength(ArPruJ,Gn); SetLength(ArPruI,Gn);SetLength(ArPruN,Gn); SetLength(RozptylJ,Gn); SetLength(RozptylI,Gn); SetLength(SmerOdI,Gn); SetLength(SmerOdJ,Gn);
  ProgressBar1.Min := 0;
  ProgressBar1.Max :=Gn ;
  ProgressBar1.Step := 1;//round(ProgressBar1.width/height);
  ProgressBar1.Position := 0;
  for b := 0 to Gn-1 do
  begin
    ProgressBar1.StepIt;
    ProgressBar1.Update;
    Data:=TList.Create;
    MaxJas[b]:=0; ArPruI[b]:=0; ArpruJ[b]:=0;ArpruN[b]:=0; RozptylI[b]:=0; RozptylJ[b]:=0;
    for I := 0 to GLPixelData[b].Count-1 do  //Najdeme maximální jas (normu vektoru posuvu) v b. obrázku
    begin
      PPixel:=GLPixelData[b].Items[I];
      ArPruI[b]:= ArPruI[b]+PPixel.PosunI;
      ArpruJ[b]:=ArPruJ[b]+PPixel.PosunJ;
      ArPruN[b]:=ArPruN[b]+sqrt(sqr(PPixel.PosunI)+sqr(PPixel.PosunJ));
      Jas:=sqrt(PPixel.PosunI*PPixel.PosunI+PPixel.PosunJ*PPixel.PosunJ);
      if MaxJas[b]<Jas then MaxJas[b]:=Jas;
    end;
    ArPruI[b]:=ArPruI[b]/(GLPixelData[b].Count);
    ArPruJ[b]:=ArPruJ[b]/(GLPixelData[b].Count);
    ArPruN[b]:=ArPruN[b]/(GLPixelData[b].Count);
    for I := 0 to GLPixelData[b].Count-1 do
    begin
      PPixel:=GLPixelData[b].Items[I];
      RozptylJ[b]:=RozptylJ[b]+sqr(PPixel.PosunJ-ArPruJ[b]);
      RozptylI[b]:=RozptylI[b]+sqr(PPixel.PosunI-ArPruI[b]);
      ar:=arg(PPixel.PosunJ,PPixel.PosunI);   //Posun J = posun ve smìru reálné osy, Posun I = posun ve smìru imaginární osy, spoète úhel od reálné osy
      Data.Add(MoveVector(ar,MaxJas[b],sqrt(PPixel.PosunI*PPixel.PosunI+PPixel.PosunJ*PPixel.PosunJ)));  //Nahraje RGB barvy, jimi si má vykreslit pøíslušnı pixel do listu dat
    end;
    RozptylJ[b]:=RozptylJ[b]/(GLPixelData[b].Count);
    RozptylI[b]:=RozptylI[b]/(GLPixelData[b].Count);
    SmerOdJ[b]:=sqrt(RozptylJ[b]);
    SmerOdI[b]:=sqrt(RozptylI[b]);
    str:=GArrayOfNames[b];//extractfilename(FileListBox1.Items[b]);  //extrahuje název
    SetLength(str, Length(str) - 4);   //odstraní koncovku .tif
    BitmapPosuvy(GWidth,GHeight,Data,'MapaPosuvu'+str);//Vykreslí a uloí posuvy
    for I := 0 to GLPixelData[0].Count-1 do dispose(Data.Items[I]); //smae pamì vyhrazenou pro posuvy
    Data.free;
    WriteLn(statFile,GArrayOfNames[b]+';'+floattostr(ArPruJ[b])+';'+floattostr(ArPruI[b])+';'+floattostr(ArPruN[b])+';'+floattostr(RozptylJ[b])+';'+floattostr(RozptylI[b])+';'+floattostr(SmerOdJ[b])+';'+floattostr(SmerOdI[b])+';'+floattostr(MaxJas[b]));
  end;
  CloseFile(statFile);
  //Vypíše statistická data do stringridu
  StringGrid1.ColCount:=9;StringGrid1.RowCount:=Gn+1;
  StringGrid1.Cells[1,0]:='Ar. prùmìr posuvù ve vodorovném smìru';StringGrid1.Cells[2,0]:='Ar. prùmìr posuvù ve svislém smìru'; StringGrid1.Cells[3,0]:='Ar. prùmìr normy posuvù';
  StringGrid1.Cells[4,0]:='Rozptyl posuvù ve vodorovném smìru';StringGrid1.Cells[5,0]:='Rozptyl posuvù ve svislém smìru';
  StringGrid1.Cells[6,0]:='Smìrodatná odchylka posuvù ve vodorovném smìru';StringGrid1.Cells[7,0]:='Smìrodatná odchylka posuvù ve svislém smìru';
  StringGrid1.Cells[8,0]:='Nejvìtší posuv ve smyslu normy';
  for I := 0 to Gn-1  do
  begin
    StringGrid1.Cells[0,I+1]:=GArrayOfNames[I]; StringGrid1.Cells[1,I+1]:=floattostr(ArPruJ[I]);
    StringGrid1.Cells[2,I+1]:=floattostr(ArPruI[I]); StringGrid1.Cells[3,I+1]:=floattostr(ArPruN[I]);
    StringGrid1.Cells[4,I+1]:=floattostr(RozptylJ[I]);
    StringGrid1.Cells[5,I+1]:=floattostr(RozptylI[I]); StringGrid1.Cells[6,I+1]:= floattostr(SmerOdJ[I]);
    StringGrid1.Cells[7,I+1]:=floattostr(SmerOdI[I]); StringGrid1.Cells[8,I+1]:=floattostr(MaxJas[I]);
  end;
  SetLength(MaxJas,0); SetLength(ArPruI,0); SetLength(ArPruJ,0);SetLength(ArPruN,0);SetLength(RozptylI,0);SetLength(RozptylJ,0);SetLength(SmerOdJ,0);SetLength(SmerOdI,0);
  showmessage('Posuvy vykresleny.');
//  ReportMemoryLeaksOnShutdown := True;
end;


procedure TMainForm.BilinearniInterpolaceClick(Sender: TObject);
var I,J,width,height,size,max,b,K,prumer,Sum:integer;
    Data:PPoleWord;
    B2: array of PPoleWord;
    PPixel:GPPointsRecordKor;
    str:string;
    MaxJasPrumer,MaxJasHrany:single;
    GPruObrazOstr: PPoleWord;
    M:array of array of ShortInt;
    dimM:byte;
    r,s:ShortInt;
    Grad:PSingleArray;
begin
  BilinearniInterpolace.enabled:=false;
  CreateDir(ExtractFilePath(ParamStr(0))+'\'+GNameOfFolder+'\'+'UpravenáData');
  width:=Gwidth; height:=Gheight; size:=GSize;
  SetLength(B2,Gn);
  for b := 0 to Gn-1 do
  begin
    GetMem(B2[b],Size);GetMem(Data,Size);max:=0;
    for I:=0 to (height-1) do
    begin
      for J:=0 to (width-1) do
      begin
        PPixel:=GLPixeldata[b].Items[I*width+J];
        Data^[I*width+J]:=BiInterpolation(J+PPixel^.PosunJ,I+PPixel^.PosunI,width,b);
        if (Data^[I*width+J]>max) then max:=Data^[I*width+J];
      end;
    end;
    str:=GArrayOfNames[b];//extractfilename(FileListBox1.Items[b]);
    SetLength(str, Length(str) - 4);   //odstraní koncovku .tif
    NewBitmapUprData(width,height,max,Data,'UpravenyObraz'+str);
    Move(Data^,B2[b]^,Size);//Zkopíruje originální data do nového pole B2[q]^
    for I:=0 to (GLPixelData[b].Count-1) do Dispose(GLPixelData[b].Items[I]);   //smae jdnotlivé pointery
    GLPixelData[b].Free;  //smae list vrcholù trojúhelníkové sítì
    FreeMem(Data);//Uvolní pointer Data
  end;
  SetLength(GLPixelData,0);
  MaxJasPrumer:=0;                      //MaxJ = Maximální jas pixelu ve zprùmìrovaném obraze hran
  for I:=0 to (height-1) do
  begin
    for J:=0 to (width-1) do
    begin
      prumer:=0;
      for K := 0 to (Gn-1) do
        begin
          prumer:=prumer+B2[K]^[I*width+J];
        end;
      GPruObraz[I*width+J]:=round((prumer/Gn)); // prùmìrnı jas pixelu [I*width+J] ve zprùmìrovaném obraze
      if(GPruObraz[I*width+J]>MaxJasPrumer) then MaxJasPrumer:=GPruObraz[I*width+J]; //Zjistí maximální jas pixelu ve zprùmìrovaném obraze
    end;
  end;
  NewBitmap(width,height,MaxJasPrumer,GPruObraz,'ZprumerovanyObrazUpraveny');
  //ZAOSTØENÍ ZPRÙMÌROVANÉHO OBRAZU
  GetMem(GPruObrazOstr,GSize);
  dimM:=3;
  SetLength(M, dimM);         //nastaví poèet øádkù masky (matice)
  for I := 0 to Length(M)-1 do SetLength(M[I], dimM);//kadı øádek má stejnı poèet sloupcù, jako má matice M poèet øádkù (vytvoøíme ètvercovou matici)
  M[0,0]:=-1; M[0,1]:=-1; M[0,2]:=-1;
  M[1,0]:=-1; M[1,1]:=9; M[1,2]:=-1;
  M[2,0]:=-1; M[2,1]:=-1; M[2,2]:=-1;
  for I := 1 to GHeight-2 do
  begin
    for J := 1 to GWidth-2 do
    begin
      Sum:=0;
      for r := -1 to dimM-2 do
      begin
        for s := -1 to dimM-2 do
        begin
          Sum:=Sum+M[r+1,s+1]*GPruObraz[(I+r)*width+(J+s)];
        end;
      end;
      GPruObrazOstr[I*width+J]:=Sum;
    end;
  end;
  MaxJasPrumer:=0;                      //MaxJ = Maximální jas pixelu ve zprùmìrovaném obraze hran
  for I:=0 to (height-1) do
  begin
    for J:=0 to (width-1) do
    begin
      if (I=0) or (I=height-1) or (J=0) or (J=width-1) then GPruObrazOstr[I*width+J]:=GPruObraz[I*width+J];
      if(GPruObrazOstr[I*width+J]>MaxJasPrumer) then MaxJasPrumer:=GPruObrazOstr[I*width+J]; //Zjistí maximální jas pixelu ve zprùmìrovaném obraze
    end;
  end;
  NewBitmap(width,height,MaxJasPrumer,GPruObrazOstr,'ZprumerovanyObrazUpravenyZaostreny');
  //Vytvoøení obrazu hran ze zprumerovaneho upraveneho obrazu
  GetMem(Grad, width*height*4);
  MaxJasHrany:=Gradient(height,width,Size,GPruObrazOstr,Grad); // Ze zprùmìrovaného obrazu v GPruObraz se napoèítají hrany a ty jsou pak uloeny do Grad
  NewBitmap(width,height,MaxJasHrany,Grad,'Hrany_ZprumerovanyObrazUpravenyZaostreny');   //Uloí data v Grad jako bmp
  FreeMem(Grad);
   //Uvolnìní pamìti
  for I := 0 to Length(M)-1 do SetLength(M[I], 0);
  SetLength(M,0);
  FreeMem(GPruObraz,GSize); FreeMem(GPruObrazOstr);
  for I := 0 to Gn-1 do
  begin
    FreeMem(B2[I]);FreeMem(GB2[I],GSize);
  end;
  SetLength(B2,0);Setlength(GB2,0);
  GArrayOfNames:=nil; Setlength(GArrayOfNames,0);
  stringgrid1.RowCount:=0; Stringgrid1.ColCount:=0;
  Stringgrid1.Free;
  Showmessage('Bilineární interpolace dokonèena.');
  //ReportMemoryLeaksOnShutdown := True;
  //restartuje program
  ShellExecute(Handle, nil, PChar(Application.ExeName), nil, nil, SW_SHOWNORMAL);
  Close;
end;

procedure TMainForm.MaxPolChange(Sender: TObject);
begin
CheckEntry(MaxPol);
end;

procedure TMainForm.MeritkoMapyPosuvuClick(Sender: TObject);
var arg,x:integer;
    Pixels:PRGBTripleArray;
    Bitmap:TBitmap;
    Popis:textFile;
begin
  MeritkoMapyPosuvu.enabled:=false; BilinearniInterpolace.Enabled:=true;
  AssignFile(Popis,ExtractFilePath(ParamStr(0))+GNameOfFolder+'\'+'MapyPosuvù'+'\'+'PopisMìøítkaPosuvù.txt');
  ReWrite(Popis);
  WriteLn(Popis,'                 Popsání mìøítka');
  WriteLn(Popis,'Poèátek je umístìn v levém dolním rohu obrázku.');
  WriteLn(Popis,'Vodorovná osa je orientovaná zleva doprava, pøièem maximální hodnota na této ose je rovna 255.');
  WriteLn(Popis,'Hodnoty na této ose znázoròují jas pixelu. Platí, e nejjasnìjší pixely jsou ty, které dosahují maximální hodnotu posunutí ve smyslu normy vektoru.');
  WriteLn(Popis,'Svislá osa je orientována od poèátku smìrem vzhùru, pøièem maximální hodnota je rovna 359.');
  WriteLn(Popis,'Hodnoty na této ose znázoròují úhel posunutí daného pixelu.');
  WriteLn(Popis,'Pomocí hodnoty na svislé ose a hodnoty na vodorovné ose, tedy lze urèit posunutí daného pixelu.');
  WriteLn(Popis,'Maximální hodnotu posuvu nalezneme v souboru StatistickéVyhodnocení.txt.');
  CloseFile(Popis);
  Bitmap:=CreateBmp(256,360,24);
  for arg := 0 to 359 do
  begin
    Pixels:=Bitmap.ScanLine[359-arg];
    for x := 0 to 255 do
    begin
      //R sloka
      if (arg<=60) or (arg>=300) then Pixels[x].rgbtRed:=round(x)
      else if (arg>60) and (arg<120) then Pixels[x].rgbtRed:=round(-(x/60)*arg+2*x)
      else if (arg>240) and (arg<300) then Pixels[x].rgbtRed:=round((x/60)*arg-4*x)
      else Pixels[x].rgbtRed:=0;
      //G sloka
      if (arg>=0) and (arg<60) then Pixels[x].rgbtGreen:=round((x/60)*arg)
      else if (arg>180) and (arg<240) then Pixels[x].rgbtGreen:=round(-(x/60)*arg+4*x)
      else if (arg>=60) and (arg<=180) then Pixels[x].rgbtGreen:=round(x)
      else Pixels[x].rgbtGreen:=0;
      //B sloka
      if (arg>120) and (arg<180) then Pixels[x].rgbtBlue:=round((x/60)*arg-2*x)
      else if (arg>300) and (arg<=359) then Pixels[x].rgbtBlue:=round(-(x/60)*arg+6*x)
      else if (arg>=180) and (arg<=300) then Pixels[x].rgbtBlue:=round(x)
      else Pixels[x].rgbtBlue:=0;
    end;
  end;
  Bitmap.SaveToFile(ExtractFilePath(ParamStr(0))+GNameOfFolder+'\'+'MapyPosuvù'+'\'+'MeritkoPosuvu'+'.bmp');
  Bitmap.Free;
  showmessage('Vykreslení dokonèeno.');
 // ReportMemoryLeaksOnShutdown := True;
end;

procedure TMainForm.UlozPosuvyClick(Sender: TObject);
var I,J,b:Integer;
    txtFile:TextFile;
    Name,Posuv,Posun:string;
    PPixel:GPPointsRecordKor;
begin
  UlozPosuvy.Enabled:=false;
  //Vypíše posuvy do txt souboru
  AssignFile(txtFile,ExtractFilePath(ParamStr(0))+GNameOfFolder+'\'+'Posuv.txt');
  ReWrite(txtFile);
  name:='';Posuv:='';
  for I := 0 to Gn-1 do
  begin
    name:=name+GArrayOfNames[I]+';;';
    Posuv:=Posuv+'Vodorovnı posuv; Svislı posuv;';
  end;
  WriteLn(txtFile,Name);
  WriteLn(txtFile,Posuv);
  ProgressBar1.Position := 0;
  ProgressBar1.Update;
  ProgressBar1.Min := 0;
  ProgressBar1.Step := 1;
  ProgressBar1.Max :=Gheight ;
  for I := 0 to GHeight-1 do
  begin
    ProgressBar1.StepIt;
    ProgressBar1.Update;
    for J := 0 to GWidth-1 do
    begin
      Posun:='';
      for b := 0 to Gn-1 do
      begin
       PPixel:=GLPixelData[b].Items[I*Gwidth+J];
       Posun:=Posun+FloatToStr(PPixel^.PosunJ)+';'+FloatToStr(PPixel^.PosunI)+';';
      end;
      WriteLn(txtFile,Posun);
    end;
  end;
  CloseFile(txtFile);
  showmessage('Data byla uloena.');
  //ReportMemoryLeaksOnShutdown := True;
end;

procedure TMainForm.KrokChange(Sender: TObject);
var S:integer;
begin
  if (TrystrToInt(Krok.Text,S)) and (S>0) and (Krok.Text[1]<>'0') then
  begin
    if 5*strtoint(Krok.Text)<25 then MaxPol.Text:='25'
      else MaxPol.Text:=IntToStr(3*strtoint(Krok.Text));
    Krok.Color:=clGreen; Krok.Font.Color:=clWhite;
  end
  else
  begin
    ShowMessage('Špatnı vstup. Vstupení parametr musí bıt pøirozené èíslo');
    Krok.Color:=clRed;
  end;
end;

procedure TMainForm.KrokPomSitChange(Sender: TObject);
begin
CheckEntry(KrokPomSit);
end;

procedure TMainForm.PolomerChange(Sender: TObject);
begin
  CheckEntry(Polomer);
end;

procedure TMainForm.PolomerKorChange(Sender: TObject);
begin
CheckEntry(PolomerKor);
end;

procedure TMainForm.PolomerPomSitChange(Sender: TObject);
begin
CheckEntry(PolomerPomSit);
end;

procedure TMainForm.PrahHodChange(Sender: TObject);
begin
  CheckEntry(PrahHod);
end;

procedure TMainForm.ResetujClick(Sender: TObject);
begin
ShellExecute(Handle, nil, PChar(Application.ExeName), nil, nil, SW_SHOWNORMAL);
  Close;
end;

procedure TMainForm.OkoliKorelaceChange(Sender: TObject);
begin
CheckEntry(OkoliKorelace);
end;



end.
