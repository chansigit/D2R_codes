import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;


public class D2RStream {
    private ArrayList<String> seqs;
    private SeqInfo seqInfo;
    public D2RStream(SeqInfo i) {
        seqInfo= i;
        seqs = new ArrayList<>();
    }

    public void add(String seq) {
        seqs.add(seq);
    }


    private Float prob(String w){
        Float pr=1.0f;
        for (int i=0;i!=w.length();++i) {
            char ch = w.charAt(i);
            if (ch=='T'||ch=='C'||ch=='G'||ch=='A'||
                    ch=='t'||ch=='c'||ch=='g'||ch=='a'){
                pr *= seqInfo.letterPr[(int) Character.toLowerCase(ch)];
            }else{
                pr *=1.0f;
            }

        }
        return pr;
    }

    private HashMap<String, Integer> getKmers(String seq, int k) {
        HashMap<String, Integer> hMap = new HashMap<>();

        for (int i = 0; i < seq.length() - k + 1; ++i) {
            String kmer = seq.substring(i, i + k);
            if (!hMap.containsKey(kmer)) {
                hMap.put(kmer, 1);
            } else {
                hMap.put(kmer, hMap.get(kmer) + 1);
            }
        }
        return hMap;
    }


    public Float D2R(String seq, int k) {
        float n = seq.length() - k + 1;
        float result = 0.0f;
        HashMap<String, Integer> hMap = getKmers(seq, k);

        for (Map.Entry<String, Integer> entry : hMap.entrySet()) {
            String kmer = entry.getKey();
            float cnt = (float) entry.getValue();
            float kmerProb = prob(kmer);
            float e = n*n*kmerProb*kmerProb;
            result += cnt*(cnt-1) - e;
        }
        //System.out.println(result);
        return result/(n*(n-1));
    }

    public ArrayList<Float> runD2R(int kmerLen) {
        return seqs.parallelStream()
                .map(s -> (float) (D2R(s, kmerLen)))
                .collect(Collectors.toCollection(ArrayList::new));
    }

}


