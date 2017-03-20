import java.io.*;
import java.util.*;

//Usage: $java starving nmr_of_chm nmr_of_spc nmr_of_steps nmr_of_environments nmr_of_reactions param.csv param_reaction.csv param_metabolism.csv param_chm.csv filename_suffix step_prefix

class starving{
	public static void main(String[] args){
		int nmr_of_chm = new Integer(args[0]).intValue();
		int nmr_of_spc = new Integer(args[1]).intValue();
		int nmr_of_reactions = new Integer(args[2]).intValue();
		int nmr_of_environments = new Integer(args[3]).intValue();
		int nmr_of_steps = new Integer(args[4]).intValue();
		File param_file = new File(args[5]);
		File reaction_file = new File(args[6]);
		File metabolism_file = new File(args[7]);
		File chm_file = new File(args[8]);
		String fix = args[9];
		String step_fix = args[11] + "_";
		Boolean stable;
		if (args[10].equals("stable")){
			stable = true;
		}else{
			stable = false;
		}
		
		Double[][] params = double_table(param_file, nmr_of_spc, 5);
		System.out.println("param file loaded");
		Integer[][] reactions = int_table(reaction_file, nmr_of_reactions,2);
		System.out.println("reaction file loaded");
		Double[][] metabolisms = double_table(metabolism_file, nmr_of_spc, nmr_of_reactions);
		System.out.println("metabolism file loaded");
		Double[][] _chemical = double_table(chm_file, nmr_of_chm, 1);
		System.out.println("chm file loaded");
		Double[] chemical = new Double[nmr_of_chm];
		for (int i=0;i<nmr_of_chm;i++){
			chemical[i] = _chemical[i][0];
		}

		Double[] abundance_dist = new Double[nmr_of_spc];
		for (int i=0;i<nmr_of_spc;i++){
			abundance_dist[i] = params[i][0];
		}

		Environment[] env = new Environment[nmr_of_environments];

		Double[][][] pathways = make_pathways(reactions, metabolisms, nmr_of_spc, nmr_of_chm, nmr_of_reactions);

		/*if (stable){
			Double[] zeros = new Double[nmr_of_chm];
			for (int c=0;c<nmr_of_chm;c++){
				zeros[c] = 0.0;
			}
			for (int e=0;e<nmr_of_environments;e++){
				env[e] = new Environment(zeros);
			}
		}else{*/
			for (int e=0;e<nmr_of_environments;e++){
				env[e] = new Environment(chemical);
			}
		//}
		for (int e=0;e<nmr_of_environments;e++){
			env[e].set_spcies(nmr_of_spc, abundance_dist, params, pathways);
		}
		System.out.println("species set");

		System.out.print("running...");
		run(nmr_of_spc, nmr_of_chm, nmr_of_environments, nmr_of_steps, env, pathways, fix, step_fix, stable, chemical, abundance_dist);
	}

	private static Double[][][] make_pathways(
		Integer[][] reactions, Double[][] metabolisms, int nmr_of_spc, int nmr_of_chm, int nmr_of_reactions){
		Double[][][] out = new Double[nmr_of_spc][nmr_of_chm][nmr_of_chm];
		for (int s=0;s<nmr_of_spc;s++){
			Double[][] pathway = new Double[nmr_of_chm][nmr_of_chm];
			for (int i=0;i<nmr_of_chm;i++){
				for (int j=0;j<nmr_of_chm;j++){
					pathway[i][j] = 0.;
				}
			}
			Double[] metabolism = metabolisms[s];
			for (int m=0;m<nmr_of_reactions;m++){
				if (metabolism[m] > 0.){
					int input_index = reactions[m][0];
					int output_index = reactions[m][1];
					pathway[input_index][output_index] += metabolism[m];
				}
			}
			for (int i=0;i<pathway.length;i++){
				pathway[i] = revise(pathway[i]);
			}
			//pathway = pathway_return_chemical(pathway);
			out[s] = pathway;
		}
		return out;
	}

	private static Double[] revise(Double[] array){
		Double[] out = new Double[array.length];
		double sum = double_sum(array);
		if (sum == 0.){
			out = array;
		}else{
			for (int i=0;i<array.length;i++){
				out[i] = array[i] / sum;
			}
		}
		return out;
	}

	private static Double[][] pathway_return_chemical(Double[][] pathway){
		Double[][] out = pathway;
		for (int i=0;i<pathway.length;i++){
			for (int j=0;j<pathway.length;j++){
				out[i][j] = out[i][j] * 0.5;
			}
		}		
		for (int i=0;i<pathway.length;i++){
			out[i][i] += 0.5;
		}
		return out;
	}

	private static double double_sum(Double[] array){
		double out = 0.;
		for (int i=0;i<array.length;i++){
			out += array[i];
		}
		return out;
	}

	private static void run(int nmr_of_spc, int nmr_of_chm, int nmr_of_environments, int nmr_of_steps, Environment[] env, Double[][][] pathways, String fix, String step_fix, boolean stable, Double[] chemical, Double[] abundance_dist){
		Double[] a_chemical = new Double[nmr_of_chm];
		for (int c=0;c<nmr_of_chm;c++){
			a_chemical[c] = chemical[c];
		}
		try{
			File file = new File("result_" + fix + ".csv");
			File file2 = new File("result_chm_" + fix + ".csv");
			PrintWriter output = new PrintWriter(new BufferedWriter(new FileWriter(file)));
			PrintWriter output2 = new PrintWriter(new BufferedWriter(new FileWriter(file2)));
			String head = "spc";
			String head2 = "chm";
			for(int s=0;s<nmr_of_spc;s++){
				head += ",spc"+s+"_abundance_total";
			}
			for (int c=0;c<nmr_of_chm;c++){
				head2 += ",chm"+c+"_abundance";
			}
			output.println(head);
			output2.println(head2);

			for (int i=0;i<nmr_of_steps;i++){
				for (int e=0;e<nmr_of_environments;e++){
					env[e].take_step(pathways);
					if (stable){
						env[e].add_chemical(a_chemical, 0.001);
					}
					if (i==100000){
						for (int c=30;c<50;c++){
						env[e].chemical[c] += 1.;}
						//env[e].add_spc(abundance_dist);
						env[e].kill_spc(0.001);}

				}
				if (i % 100 == 99){
					int step = i+1;
					System.out.print("\rrunning... " + step + " / " + nmr_of_steps + " steps done");
					Double[][] each_abundance = new Double[nmr_of_environments][nmr_of_spc];
					Double[][] chm_abundance = new Double[nmr_of_environments][nmr_of_chm];
					for (int e=0;e<nmr_of_environments;e++){
						each_abundance[e] = env[e].abundances();
						chm_abundance[e] = env[e].chemical;
					}
					Double[] total_abundance = mean_abundance(each_abundance);
					Double[] total_chm = mean_abundance(chm_abundance);
					output.print(step_fix + i);
					output2.print(step_fix + i);
					for(int s=0;s<nmr_of_spc;s++){
						output.print("," + total_abundance[s]);
					}
					for(int c=0;c<nmr_of_chm;c++){
						output2.print("," + total_chm[c]);
					}
					output.print("\n");
					output2.print("\n");
				}
			}
			output.close();
			output2.close();
			System.out.println("\nsimulation done.");

		}catch(FileNotFoundException e){
			System.out.println(e);
		}catch(IOException e){
			System.out.println(e);
		}
	}

	private static Double[] mean_abundance(Double[][] abundances){
		Double[] out = new Double[abundances[0].length];
		double r = 1. / (double)abundances.length;
		for (int i=0;i<out.length;i++){
			out[i] = 0.;
		}
		for (int i=0;i<abundances.length;i++){
			for (int j=0;j<out.length;j++){
				out[j] += abundances[i][j] * r;
			}
		}
		return out;
	}

	private static Double[][] double_table(File file, int nmr_of_indices, int nmr_of_columns){
		Double[][] out = new Double[nmr_of_indices][nmr_of_columns];
		try{
			BufferedReader br = new BufferedReader(new FileReader(file));
			String str;
			int l = -1;
			while((str = br.readLine()) != null){
				if (l > -1){
					String[] values = str.split(",", 0);
					for (int i=0;i<nmr_of_columns;i++){
					 	out[l][i] = new Double(values[i+1]).doubleValue();
					}
				}
				l++;
			}
			br.close();
		}catch(FileNotFoundException e){
			System.out.println(e);
		}catch(IOException e){
			System.out.println(e);
		}
		return out;
	}

	private static Integer[][] int_table(File file, int nmr_of_indices, int nmr_of_columns){
		Integer[][] out = new Integer[nmr_of_indices][nmr_of_columns];
		try{
			BufferedReader br = new BufferedReader(new FileReader(file));
			String str;
			int l = -1;
			while((str = br.readLine()) != null){
				if (l > -1){
					String[] values = str.split(",", 0);
					for (int i=0;i<nmr_of_columns;i++){
					 	out[l][i] = new Integer(values[i+1]).intValue();
					}
				}
				l++;
			}
			br.close();
		}catch(FileNotFoundException e){
			System.out.println(e);
		}catch(IOException e){
			System.out.println(e);
		}
		return out;
	}

}


class Environment{
	public static Gamma gmm = new Gamma();
	public Double[] chemical;
	public Species[] spc;
	private double total_abundance;
	private static final double clock = 0.01;
	private static final double hill_trans = 1.;
	private static final double hill_growth = 1.;
	private static final double hill_death = 2.;
	private static final double death_k = 0.01;
	private static double chm_spc_rate;
	private static final double capacity = 1.;
	private static final double abundance_limit = 0.0000000000001;
	private static final double sampling_size = 10.; ///changed
	

	Environment(Double[] set_chemical){
		chemical = new Double[set_chemical.length];
		for (int c=0;c<set_chemical.length;c++){
			chemical[c] = set_chemical[c];
		}
		//chm_spc_rate = (double)set_chemical.length;
		chm_spc_rate = 1.;
	}

	public void set_spcies(int nmr_of_spc, Double[] abundance_dist, Double[][] params, Double[][][] pathways){
		spc = new Species[nmr_of_spc];
		double __total_abundance = double_sum(abundance_dist);
		for (int s=0;s<spc.length;s++){
			ArrayList<Integer> set_mii = make_metabolism_input_indices(pathways[s]);
			double sample_abundance = abundance_dist[s] / __total_abundance * sampling_size;
			double sigma = Math.sqrt(1.0 / sample_abundance);
			//spc[s] = new Species(sample(abundance_dist[s], sigma), params[s][1], params[s][2], params[s][3], params[s][4], set_mii);
			spc[s] = new Species(abundance_dist[s], params[s][1], params[s][2], params[s][3], params[s][4], set_mii);
		}
		total_abundance = double_sum(abundances());
	}

	public ArrayList<Integer> make_metabolism_input_indices(Double[][] pathway){
		ArrayList<Integer> out = new ArrayList<Integer>();
		for (int c=0;c<pathway.length;c++){
			if ( double_sum(pathway[c]) > 0.){
				out.add(c);
			}
		}
		return out;
	}

	private static double double_sum(Double[] array){
		double out = 0.;
		for (int i=0;i<array.length;i++){
			out += array[i];
		}
		return out;
	}

	private static double sample(double value, double sigma){
		double kappa = 1. / sigma;
		double theta = sigma;
		// mean = 1.0, variance = sigma
		double out = value * gmm.random_Gamma(kappa,theta);
		return out;
	}

	public void take_step(Double[][][] pathways){
		List<Integer> index = new ArrayList<>();
		for (int i=0;i<spc.length;i++){index.add(i);}
		Collections.shuffle(index);
		for(int s=0;s<spc.length;s++){
			if (spc[index.get(s)].alive){
				move_one_species(index.get(s), pathways);
				//__move_one_species(index.get(s), pathways);
			}
		}
	}

	private void move_one_species(int spc_index, Double[][][] pathways){
		/*
		if (total_abundance > capacity || total_abundance < 0.){
			System.out.println("\n"+total_abundance);
			System.exit(1);
		}
		*/
		for (int array_index=0;array_index<spc[spc_index].metabolism_input_indices.size();array_index++){
			int chm_index = spc[spc_index].metabolism_input_indices.get(array_index);
			Double[] pathway = pathways[spc_index][chm_index];
			move_one_active(spc_index, chm_index, array_index, pathway);
		}
		int chm_index = max_chm_index(spc[spc_index]);
		int array_index = spc[spc_index].metabolism_input_indices.indexOf(chm_index);
		Double[] pathway = pathways[spc_index][chm_index];
		move_inactive_cells(spc_index, chm_index, array_index, pathway);
		extinction(spc[spc_index]);
	}

	private void move_one_active(int spc_index, int chm_index, int array_index, Double[] pathway){
		double gr = monod(chemical[chm_index], spc[spc_index].growth_k);
		double tr = monod(chemical[chm_index], spc[spc_index].trans_k);
		//double tr = hill(spc[spc_index].active_abundance.get(array_index), 1.0 , hill_trans);
		//double ga = gr * spc[spc_index].sa * spc[spc_index].active_abundance.get(array_index) * (1. - total_abundance / capacity) * clock;
		double ga = gr * spc[spc_index].sa * spc[spc_index].active_abundance.get(array_index) * clock;
		double ta = (1. - tr) * spc[spc_index].sa * spc[spc_index].active_abundance.get(array_index) * clock;

		chemical[chm_index] -= ga * chm_spc_rate;
		if (chemical[chm_index] < 0.){
			ga += chemical[chm_index] / chm_spc_rate;
			chemical[chm_index] = 0.;
		}

		double gain_a = ga - ta;
		double abundance_a = spc[spc_index].active_abundance.get(array_index) + gain_a;

		double gain_i;
		if (abundance_a < 0.){
			gain_i = abundance_a - gain_a;
			abundance_a = 0.;
		}else{
			gain_i = ta;
		}
		
		spc[spc_index].active_abundance.set(array_index, abundance_a);
		spc[spc_index].inactive_abundance += gain_i;
		total_abundance += gain_a + gain_i;
		if (total_abundance < 0.){System.out.println("active");}

		for(int c=0;c<chemical.length;c++){
			chemical[c] += ga * pathway[c] * chm_spc_rate;
		}

	}

	private void move_inactive_cells(int spc_index, int chm_index, int array_index, Double[] pathway){
		double gr = monod(chemical[chm_index], spc[spc_index].growth_k);
		//double dr = 1. / sampling_size ;
		double dr = 1. - gr;
		//double dr = hill(spc[spc_index].inactive_abundance, death_k, hill_death);
		//double gi = gr * spc[spc_index].si * spc[spc_index].inactive_abundance * (1. - total_abundance / capacity) * clock;
		double gi = gr * spc[spc_index].si * spc[spc_index].inactive_abundance * clock;
		double di = dr * spc[spc_index].si * spc[spc_index].inactive_abundance * clock;
		
		chemical[chm_index] -= gi * chm_spc_rate;
		if (chemical[chm_index] < 0.){
			gi += chemical[chm_index] / chm_spc_rate;
			chemical[chm_index] = 0.;
		} 

		double gain_a = gi;

		double abundance_a = spc[spc_index].active_abundance.get(array_index) + gain_a;
		spc[spc_index].active_abundance.set(array_index, abundance_a);
		spc[spc_index].inactive_abundance -= di;
		if (spc[spc_index].inactive_abundance < 0.){
			di -= spc[spc_index].inactive_abundance;
			spc[spc_index].inactive_abundance = 0.;
		}
		total_abundance += gain_a - di;
		if (total_abundance < 0.){System.out.println("inactive");}

		for(int c=0;c<chemical.length;c++){
			chemical[c] += gi * pathway[c] * chm_spc_rate;
		}
	}

	private void __move_one_species(int spc_index, Double[][][] pathways){
		int chm_index = max_chm_index(spc[spc_index]);
		int array_index = spc[spc_index].metabolism_input_indices.indexOf(chm_index);
		Double[] pathway = pathways[spc_index][chm_index];
		__move_inactive_cells(spc_index, chm_index, array_index, pathway);
		extinction(spc[spc_index]);
	}

	private void __move_inactive_cells(int spc_index, int chm_index, int array_index, Double[] pathway){
		double gr = monod(chemical[chm_index], spc[spc_index].growth_k);
		double tr = 0.;
		//double dr = Math.pow(gr, 1.2);
		//double dr = 1. / sampling_size * 10.;
		double dr = 1. - gr;
		double gi = gr * spc[spc_index].si * spc[spc_index].inactive_abundance * clock;
		//double gi = gr * spc[spc_index].si * spc[spc_index].inactive_abundance * (1. - total_abundance / capacity) * clock;
		double di = dr * spc[spc_index].si * spc[spc_index].inactive_abundance * clock;

		chemical[chm_index] -= gi * chm_spc_rate;
		if (chemical[chm_index] < 0.){
			gi += chemical[chm_index] / chm_spc_rate;
			chemical[chm_index] = 0.;
		} 

		double gain_i = gi * (1. - tr);
		//double gain_a = gi * tr;

		//double abundance_a = spc[spc_index].active_abundance.get(array_index) + gain_a;
		//spc[spc_index].active_abundance.set(array_index, abundance_a);
		spc[spc_index].inactive_abundance += gain_i - di;
		total_abundance += gain_i - di;
		for(int c=0;c<chemical.length;c++){
			chemical[c] += gi * pathway[c] * chm_spc_rate;
		}	
	}

	private void extinction(Species spc){
		double total = 0.;
		if (spc.inactive_abundance < abundance_limit){
			spc.inactive_abundance = 0.;
		}
		total += spc.inactive_abundance;
		for (int i=0;i<spc.active_abundance.size();i++){
			if (spc.active_abundance.get(i) < abundance_limit){
				spc.active_abundance.set(i, 0.);
			}
			total += spc.active_abundance.get(i);
		}
		
		if (total < abundance_limit){
			spc.alive = false;
		}
	}

	public void add_chemical(Double[] a_chemical, double rate){
		for (int i=0;i<chemical.length;i++){
			chemical[i] += a_chemical[i] * clock * rate;
		}
	}

	public void kill_spc(double rate){
		for (int s=0;s<spc.length;s++){
			spc[s].kill(rate);
		}
	}

	public void add_spc(Double[] dist){
		double __total_abundance = double_sum(dist);
		for (int s=0;s<spc.length;s++){
			double sample_abundance = dist[s] / __total_abundance * sampling_size;
			double sigma = Math.sqrt(1.0 / sample_abundance);
			spc[s].inactive_abundance += sample(dist[s],sigma) * 10.; ///changed
		}		
	}

	public Double[] abundances(){
		Double[] out = new Double[spc.length];
		for (int i=0;i<spc.length;i++){
			out[i] = spc[i].total_abundance();
		}
		return out;
	}

 /*
	private double linear(double c, double k, double s){
		if (c < k){
			return 0.;
		}else{
			return Math.min(1, s * (c - k));
		}
	}
*/

	private double monod(double c, double k){
		double out = c / ( k + c );
		return out;
	}

	private double hill(double c, double k, double n){
		double cn = Math.pow(c, n);
		double out = cn / ( Math.pow(k,n) + cn );
		return out;
	}

	private double poisson(double c, double k){
		return 1. - Math.pow(Math.E, (-1) * c / k);
	}

	private int max_chm_index(Species spc){
		int len = spc.metabolism_input_indices.size();
		double temp = 0.;
		int out = spc.metabolism_input_indices.get(0);
		for (int i=1;i<len;i++){
			int index = spc.metabolism_input_indices.get(i);
			if (chemical[index] > temp){
				temp = chemical[index];
				out = index;
			}
		}
		return out;
	}	

}

class Species{
	public double trans_k, growth_k, si, sa;
	public ArrayList<Integer> metabolism_input_indices;
	public ArrayList<Double> active_abundance;
	public double inactive_abundance;
	public boolean alive = true;

	Species(double set_abundance, double set_tk, double set_gk, double set_sa, double set_si, ArrayList<Integer> set_mii){
		inactive_abundance = set_abundance * 1.0;
		trans_k = set_tk;
		growth_k = set_gk;
		sa = set_sa;
		si = set_si;
		metabolism_input_indices = set_mii;
		active_abundance = new ArrayList<Double>();
		for (int i=0;i<metabolism_input_indices.size();i++){
			active_abundance.add(set_abundance * 0.0 / (double)set_mii.size() );
		}
	}

	public double total_abundance(){
		double total = 0.;
		for (int i=0;i<active_abundance.size();i++){
			total += active_abundance.get(i);
		}
		total += inactive_abundance;
		return total;
	}

	public void kill(double rate){
		for (int i=0;i<active_abundance.size();i++){
			active_abundance.set(i, active_abundance.get(i) * rate);
		}
		inactive_abundance = inactive_abundance * rate;
	}
}

class Gamma{
	Gamma(){}
	public static double random_Gamma(double kappa, double theta){
		Random rnd = new Random();
		if(kappa>1){
			double c1 = kappa - 0.3333333333333333;
			double c2 = 1/Math.sqrt(9*c1);
			while(true){
				double z = rnd.nextGaussian();
				while(z*c2 <= -1){
					z = rnd.nextGaussian();
				}
				double v = Math.pow((1+z*c2),3);
				double u = Math.random();
				if(u < 1-0.0331*Math.pow(z,4)){
					return c1*v*theta;
				}else if(Math.log(u) < 0.5*z*z+c1*(1-v+Math.log(v))){
					return c1*v*theta;
				}
			}
		}else{
			double c1 = 0.07+0.75*Math.sqrt(1-kappa);
			double c2 = 1+kappa/Math.pow(Math.E, c1)/c1;
			double c3 = 1/kappa;
			while(true){
				double u1 = Math.random();
				double u2 = Math.random();
				double v = c2*u1;
				if(v<1){
					double x = c1*Math.pow(v,c3);
					if(u2<(2-x)/(2+x)){
						return x*theta;
					}else if(u2<Math.pow(Math.E, -x)){
						return x*theta;
					}
				}else{
					double x = -Math.log(c1*c3*(c2-v));
					double y = x/c1;
					if(u2*(kappa+y-kappa*y)<1){
						return x*theta;
					}else if(u2<Math.pow(y,kappa-1)){
						return x*theta;
					}
				}
			}
		}
	}
}
