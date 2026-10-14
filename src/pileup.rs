/*
* IMGT/StatAssembly
Available under EUPL license
Made by: Guilhem Zeitoun and IMGT Team
*/
use extended_htslib::bam::{
    self,
    pileup::{RustPileupConfig, RustPileups},
};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::{
    fs,
    io::{
        self,
        ErrorKind::{BrokenPipe, InvalidInput, NotFound},
        Write,
    },
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

use crate::{
    givename,
    r#struct::{Args, Filecrea, LocusInfos, Positionstrand, Species},
    submissions::getprogressbarclassic,
};
pub(crate) fn checksamtools() -> io::Result<()> {
    println!("Checking samtools is present.");
    let pipe = Stdio::piped();
    let output = Command::new("samtools")
        .arg("version")
        .stdout(Stdio::null())
        .stderr(pipe)
        .output();
    match output {
        Ok(s) if s.status.success() => Ok(()),
        Ok(s) => Err(io::Error::new(
            BrokenPipe,
            format!(
                "Samtools command failed. Error is {}",
                String::from_utf8_lossy(&s.stderr)
            ),
        )),
        Err(e) => Err(io::Error::new(
            NotFound,
            format!("Samtools not found. Error is {}", e),
        )),
    }
}
pub(crate) fn commandlinepileup(
    species: &Species,
    args: &Args,
    locus: &[LocusInfos],
) -> io::Result<()> {
    checksamtools()?;
    let path = Path::join(&args.outdir, "pileups");
    match fs::create_dir(&path) {
        Ok(_) => (),
        Err(e) if e.kind() == io::ErrorKind::AlreadyExists => (), //No alert if folder already exists
        Err(e) => return Err(e),
    };
    println!("Generating pileups...");
    let progressbar = getprogressbarclassic(locus.len().try_into().unwrap_or_default());
    locus
        .par_iter()
        .map(|loc| {
            let file = Filecrea::createfrompath(Path::join(
                path.as_path(),
                PathBuf::from(
                    givename(
                        species,
                        loc.getlocus(),
                        loc.getcontig(),
                        loc.gethaplotype().isprimary(),
                        "pileup.txt",
                        false,
                    )
                    .to_string(),
                ),
            ))
            .setfile()?;
            let pipe = Stdio::piped();
            let region = format!(
                "{}:{}-{}",
                loc.contig,
                loc.getpositionstart().getobasedpos(),
                loc.getpositionend().getobasedpos()
            );
            let text = match (
                args.assembly.as_ref(),
                args.assemblyindex.as_ref(),
                args.file.as_ref(),
                args.index.as_ref(),
            ) {
                (_, _, None, _) => return Err(io::Error::new(InvalidInput, "No BAM file given")),
                (Some(a), None, Some(b), None) => vec![
                    "-f".to_string(),
                    a.display().to_string(),
                    b.display().to_string(),
                ],
                (Some(a), Some(b), Some(c), None) => vec![
                    "-f".to_string(),
                    a.display().to_string(),
                    b.display().to_string(),
                    c.display().to_string(),
                ],
                (Some(a), Some(b), Some(c), Some(d)) => vec![
                    "-f".to_string(),
                    a.display().to_string(),
                    b.display().to_string(),
                    c.display().to_string(),
                    d.display().to_string(),
                ],
                (None, None, Some(a), None) => vec![a.display().to_string()],
                (None, None, Some(a), Some(b)) | (None, Some(_), Some(a), Some(b)) => {
                    vec![a.display().to_string(), b.display().to_string()]
                }
                (None, Some(_), Some(a), None) => vec![a.display().to_string()],
                (Some(a), None, Some(b), Some(c)) => vec![
                    "-f".to_string(),
                    a.display().to_string(),
                    b.display().to_string(),
                    c.display().to_string(),
                ],
            };
            let mut args: Vec<String> = vec![
                "mpileup",
                "--reverse-del",
                "-a",
                "-B",
                "-Q",
                "0",
                "-q",
                "0",
                "-r",
                region.as_str(),
            ]
            .into_iter()
            .map(|b| b.to_string())
            .collect();
            args.extend(text);
            let _ = Command::new("samtools")
                .args(args)
                .stdout(file)
                .stderr(pipe)
                .output()?;
            if let Ok(b) = &progressbar {
                b.set_position(b.position().saturating_add(1));
            }
            Ok(())
        })
        .collect()
}
#[allow(unused)]
pub(crate) fn printpileup(species: &Species, args: &Args, locus: &[LocusInfos]) -> io::Result<()> {
    println!("Generating pileups...");
    let progressbar = getprogressbarclassic(locus.len().try_into().unwrap_or_default());
    let path = args
        .file
        .as_ref()
        .ok_or_else(|| io::Error::new(InvalidInput, "No BAM given"))?;
    locus
        .par_iter()
        .map(|loc| {
            let mut file = Filecrea::createfrompath(Path::join(
                &args.outdir,
                PathBuf::from(
                    givename(
                        species,
                        loc.getlocus(),
                        loc.getcontig(),
                        loc.gethaplotype().isprimary(),
                        "pileup.txt",
                        false,
                    )
                    .to_string(),
                ),
            ))
            .setfile()?;
            let d = RustPileups::new(
                path,
                args.index.as_ref(),
                args.assembly.as_ref(),
                args.assemblyindex.as_ref(),
                bam::FetchDefinition::RegionString(
                    loc.contig.as_bytes(),
                    loc.getpositionstart().getzbasedpos(),
                    loc.getpositionend().getzbasedpos(),
                ),
                RustPileupConfig::default(),
            )
            .map_err(|d| io::Error::new(InvalidInput, d))?;
            file.write_all(d.to_string().as_bytes())?;
            file.flush()?;
            if let Ok(p) = &progressbar {
                p.set_length(p.length().unwrap_or_default().saturating_add(1));
            }
            Ok(())
        })
        .collect::<io::Result<()>>()?;
    if let Ok(p) = progressbar {
        p.with_finish(indicatif::ProgressFinish::AndLeave);
    }
    println!("Pileups finished");
    Ok(())
}
