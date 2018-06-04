import java.util.ArrayList;

public class SeqInfo {
    public ArrayList<Integer> percentages;
    public String alphabetStr;
    public Float[] letterPr;

    public String getAlphabetStr() {
        return alphabetStr;
    }
    public void   setAlphabetStr(String alphabetStr) {
        this.alphabetStr = alphabetStr;
    }
    public ArrayList<Integer> getPercentages() {
        return percentages;
    }
    public void setPercentages(ArrayList<Integer> percentages) {
        this.percentages = percentages;
    }

    public SeqInfo(String alphabetStr, ArrayList<Integer> percentages) {
        this.alphabetStr = alphabetStr;
        this.percentages = percentages;
        letterPr = new Float[256];

        for (int i=0; i!=alphabetStr.length(); ++i){
            char ch = alphabetStr.charAt(i);
            letterPr[(int)ch] = percentages.get(i)/100.0f;
        }
    }


    public String toString(){
        StringBuffer sb = new StringBuffer();
        for (int i=0; i!=alphabetStr.length(); ++i){
            char  ch = alphabetStr.charAt(i);
            float pr = percentages.get(i)/100.0f;
            sb.append(String.format("Pr[%c]=%.2f\t",ch, pr));
        }
        return sb.toString();
    }
}
