import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;
public class D2RReadFilter {

    public static void main(String[] args) throws FileNotFoundException {
        // Usage: D2RReadFilter    reads.fa    kmersize    threshold    pt pc pg pa
        String readFile = args[0];
        int    kmersize = Integer.parseInt( args[1]) ;
        float threshold = Float.parseFloat(args[2]);
        float pt,pc,pg,pa;
        pt = Float.parseFloat(args[3]);
        pc = Float.parseFloat(args[4]);
        pg = Float.parseFloat(args[5]);
        pa = Float.parseFloat(args[6]);

        SeqInfo sInfo = new SeqInfo("tcga", new ArrayList<>(Arrays.asList(
                (int)(pt/(pt+pc+pg+pa) * 100.0f),
                (int)(pc/(pt+pc+pg+pa) * 100.0f),
                (int)(pg/(pt+pc+pg+pa) * 100.0f),
                (int)(pa/(pt+pc+pg+pa) * 100.0f))));
        //System.out.println(sInfo.toString());
        D2RStream computingCenter = new D2RStream(sInfo);
        ArrayList<String> rawSequences = new ArrayList<>();
        ArrayList<String> identifiers =new ArrayList<>();
        //========================================================================
        //========================================================================

        long startTime,endTime;
        System.out.println("Loading sequencing reads...");
        startTime=System.currentTimeMillis();   //获取开始时间
        ///boolean first = true;
        int cnt  =0 ;
        try (Scanner sc = new Scanner(new File(readFile))) {
            while (sc.hasNextLine()) {
                String line = sc.nextLine().trim();
                if (line.charAt(0) == '>') {
                    /*
                    if (first)
                        first = false;
                    else
                        System.out.println();
                    System.out.printf("%s: ", line.substring(1));
                    */
                    ;
                    identifiers.add(line);
                } else {
                    cnt++;
                    //System.out.println(line);
                    rawSequences.add(line);
                    computingCenter.add(line);
                }
            }
        }
        endTime=System.currentTimeMillis(); //获取结束时间
        System.out.printf("%d reads loaded, taking %d ms.\n",cnt,(endTime-startTime));

        //========================================================================
        //========================================================================

        System.out.println("Computing D2R statistics...");
        startTime=System.currentTimeMillis();   //获取开始时间
        ArrayList<Float> values = computingCenter.runD2R(kmersize);
        endTime=System.currentTimeMillis(); //获取结束时间
        System.out.printf("Computation finished, taking %d ms.",(endTime-startTime));

        //========================================================================
        //========================================================================
        System.out.println("Output sequences with repeats:");
        startTime=System.currentTimeMillis();   //获取开始时间
        int repCnt = 0;

        try{
            FileWriter fw = new FileWriter(readFile+"_filtered.fasta");
            for(int i=0;i!=values.size();++i){
                if (values.get(i)> threshold){
                    ///System.out.println("D2R stat of >#"+i+" ="+values.get(i));
                    //System.out.println(identifiers.get(i));
                    //System.out.println(rawSequences.get(i));
                    fw.write(identifiers.get(i)+"\n");

                    fw.write(rawSequences.get(i)+"\n");
                    repCnt++;
                }
            }
            fw.close();
        }catch (IOException e){
            System.out.println("IO Exceptions!");
        }


        endTime=System.currentTimeMillis(); //获取结束时间
        System.out.printf("Output finished, capturing %d reads out of total %d reads, taking %d ms.", repCnt,cnt,(endTime-startTime));

    }
}
