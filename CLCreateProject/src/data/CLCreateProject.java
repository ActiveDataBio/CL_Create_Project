package data;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.Arrays;

public class CLCreateProject {

  private final String baseUrl = "http://localhost:8080";//this needs to be changed when new rest apis are deployed to the server
  private static String osName = System.getProperty("os.name");
  private static final String RSCRIPT_COMMAND =
      (osName.contains("Windows")) ? "C:\\Program Files\\R\\R-3.2.1\\bin\\Rscript.exe" : "Rscript";
  private static final String RSCRIPT_DATA_IMPORT = "rdata2dendro.R";
  private static final String RSCRIPT_RUN_BATCH = "run_batch.R";

  private String organism;
  private File destination_path;
  private File rdata;
  private File metadata;
  private File background;
  private File geneset;

  /**
   * constructor performs all operations required to create a Adbio project
   * @param args
   */
  public CLCreateProject(String[] args) {
    if(args.length < 4){
      System.out.println(usage());
      System.exit(1);
    }
    this.organism = args[0];
    this.destination_path = new File(args[1]);
    if (!this.destination_path.exists())
      this.destination_path.mkdirs();
    this.rdata = createFile(new File(this.destination_path, "data.RData"));
    this.metadata = createFile(new File(this.destination_path, "metadata.tsv"));
    this.background = createFile(new File(this.destination_path, "background.csv"));
    this.geneset = createFile(new File(this.destination_path, "genesets.RData"));
    try {
      System.out.println("Transfering files to destination folder");
      if(!transferFiles(Arrays.copyOfRange(args, 2, args.length))){
        System.out.println("An error occured while transfering files to the destination folder");
        System.exit(1);
      }
      System.out.println("Generating cluster information");
      generateClusterInfo();
      System.out.println("Getting geneset file and background file");
      httpFiles();
      System.out.println("Running batch mode project generation almost complete");
      runBatchMode();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
  /**
   * function to download the files required for a project
   */
  private void httpFiles() {
    try {
      downloadFile(this.geneset, baseUrl + "/bioviz/services/projectGenesetFile/" + organism);
      if (this.background.length() == 0) {
        String proteinIds = null;
        File matrix = new File(this.destination_path, "matrix.csv");
        if (!matrix.exists()) {
          return;
        }
        String idLine = fileToString(matrix).split("\n")[0];
        proteinIds = idLine.substring(idLine.indexOf(",") + 1).replace("\"", "");
        downloadFile(this.background,
            baseUrl + "/bioviz/services/projectBackgroundFile/" + organism, proteinIds);
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
  /**
   * this method creates physical files of what's in the resource folder
   * @param str Name of the file in resources folder 
   * @return Returns the absolute path of the file created from the resource folder
   */
  private String createTempFile(String str){
    File temp = new File("./"+str);
    try {
      if(!temp.exists())
        temp.createNewFile();
    } catch (IOException e1) {
      e1.printStackTrace();
    }
    try {
      InputStream in = this.getClass().getResourceAsStream("/resources/"+str);
      FileOutputStream out = new FileOutputStream(temp);
     
      writeFileFromStream(new BufferedInputStream(in),out);
    } catch (IOException e) {
      e.printStackTrace();
    }
    temp.deleteOnExit();
    return temp.getAbsolutePath();
  }
  /**
   * This method uses a R script to generate hierarchical clustering of the data
   * @throws Exception
   */
  private void generateClusterInfo() throws Exception {
    Process p;
    p = new ProcessBuilder(RSCRIPT_COMMAND, createTempFile(RSCRIPT_DATA_IMPORT),
        this.destination_path.getAbsolutePath(), this.rdata.getAbsolutePath()).start();
    StreamGobbler errorGobbler = new StreamGobbler(p.getErrorStream(), "ERROR");
    StreamGobbler outputGobbler = new StreamGobbler(p.getInputStream(), "OUTPUT");

    errorGobbler.start();
    outputGobbler.start();
    int exitVal = p.waitFor();
    System.out.println("Generate Cluster info ExitValue: " + exitVal);
    p.destroy();

    if (exitVal != 0) {
      throw new Exception("Failed to generate dendrograms.");
    }
  }
  /**
   * this function generates the rest of the needed files using a R srcipt
   * @throws Exception
   */
  private void runBatchMode() throws Exception {
    Process p = new ProcessBuilder(RSCRIPT_COMMAND, createTempFile(RSCRIPT_RUN_BATCH),
        this.destination_path.getAbsolutePath()).start();
    StreamGobbler errorGobbler = new StreamGobbler(p.getErrorStream(), "ERROR");
    StreamGobbler outputGobbler = new StreamGobbler(p.getInputStream(), "OUTPUT");

    errorGobbler.start();
    outputGobbler.start();

    int exitVal = p.waitFor();
    System.out.println("Run Batch Mode ExitValue: " + exitVal);
    p.destroy();

    if (exitVal != 0) {
      throw new Exception("Failed to precomputing enrichment tests in batch mode.");
    }
  }
  /**
   * this function returns to the user how to order the arguments on the command line
   * @return Returns a string for using program on command line
   */
  private String usage(){
    return "CLCreateProject <abriviated organism e.g. hsa> <distination folder> <path to RData file> <path to metadata file> <optional path to background.csv";
  }

  //helpers
  /**
   * Checks if string is null or empty
   * @param str a string
   * @return true if the string is null or empty
   */
  private boolean nullOrEmpty(String str) {
    return str == null || str.isEmpty();
  }

  /**
   * This function converts a file to a string
   * @param file The file to convert to a string
   * @return Returns the contents of the file as a string 
   */
  private String fileToString(File file) {
    if (file != null && file.exists()) {
      try {
        BufferedReader br = new BufferedReader(new FileReader(file));
        StringBuilder sb = new StringBuilder();
        String line = br.readLine();
        while (line != null) {
          sb.append(line);
          sb.append(System.lineSeparator());
          line = br.readLine();
        }
        br.close();
        return sb.toString();

      } catch (IOException e) {
        e.printStackTrace();
      }
    }
    return "";
  }
  
  /**
   * This function deletes and creates a new file
   * @param file The file to recreate
   * @return Returns the file
   */
  private File createFile(File file) {
    if (file.exists())
      file.delete();
    try {
      file.createNewFile();
    } catch (IOException e) {
      e.printStackTrace();
    }
    return file;
  }

  /**
   * This function uses arguments from user input to copy user created file to the destination folder
   * @param fileStr string array of the 2 or 3 files the are required for a project
   * @return Returns true if all files successfully transfered, false otherwise.
   */
  private boolean transferFiles(String[] fileStr) {
    try {
      Files.copy(new File(fileStr[0]).toPath(), this.rdata.toPath(),
          StandardCopyOption.REPLACE_EXISTING);
      Files.copy(new File(fileStr[1]).toPath(), this.metadata.toPath(),
          StandardCopyOption.REPLACE_EXISTING);
      if (fileStr.length > 2 && !nullOrEmpty(fileStr[2])) {
        Files.copy(new File(fileStr[2]).toPath(), this.background.toPath(),
            StandardCopyOption.REPLACE_EXISTING);
      }
    } catch (IOException e) {
      e.printStackTrace();
      return false;
    }
    return true;
  }
  /**
   * This function downloads a file using post
   * @param file This is the file object the downloaded file will saved in
   * @param urlStr This is the string of the url to download the file from
   * @param proteinIds a comma separated variables string of all the protein ids needed to created a proper background.csv file
   * @throws IOException
   */
  private void downloadFile(final File file, final String urlStr, final String proteinIds)
      throws IOException {
    if (proteinIds == null)
      return;
    // Construct data
    String data =
        URLEncoder.encode("matrix", "UTF-8") + "=" + URLEncoder.encode(proteinIds, "UTF-8");
    // Send data
    URL url = new URL(urlStr);
    URLConnection conn = url.openConnection();
    conn.setDoOutput(true);
    OutputStreamWriter wr = new OutputStreamWriter(conn.getOutputStream());
    wr.write(data);
    wr.flush();
    wr.close();
    writeFileFromStream(new BufferedInputStream(conn.getInputStream()),new FileOutputStream(file));
  }
  
  /**
   * This function will download a file from a url
   * @param file The file the downloaded file will be stored in
   * @param urlString The url of the file to download
   * @throws MalformedURLException
   * @throws IOException
   */
  private void downloadFile(final File file, final String urlString)
      throws MalformedURLException, IOException {
    writeFileFromStream(new BufferedInputStream(new URL(urlString).openStream()),new FileOutputStream(file));
  }
  
  /**
   * This function copies data from an inputstream and outputs to a output stream, closes streams after
   * @param in the input stream
   * @param fout the output stream
   * @throws IOException
   */
  private void writeFileFromStream(InputStream in,OutputStream fout) throws IOException{
    try {
      final byte data[] = new byte[1024];
      int count;
      while ((count = in.read(data, 0, 1024)) != -1) {
        fout.write(data, 0, count);
      }
    } finally {
      if (in != null) {
        in.close();
      }
      if (fout != null) {
        fout.close();
      }
    }
  }

  /**
   * Main function to run on it's own
   * @param args
   */
  public static void main(String[] args) {
    System.out.println("Creating Project");
    new CLCreateProject(args);
    System.out.println("finished");
  }

}
