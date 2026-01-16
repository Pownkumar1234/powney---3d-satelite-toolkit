const fs = require('fs');
const path = require('path');

const TARGET_URL = 'https://github.com/Pownkumar1234/powney---3d-satelite-toolkit';
const SEARCH_URLS = [
    'https://www.keeptrack.space',
    'https://keeptrack.space',
    'http://keeptrack.space',
    'http://www.keeptrack.space'
];

const IGNORE_DIRS = ['node_modules', '.git', 'dist', 'build', 'coverage'];
const IGNORE_FILES = ['package-lock.json', 'yarn.lock', 'CHANGELOG.md', 'replace_links.cjs', 'analyze_links.cjs', 'analyze_links.js', 'links.json', 'found_links.txt', 'https_report.md', 'keeptrack_links.txt', 'keeptrack_links_utf8.txt'];

function walk(dir, fileList = []) {
    const files = fs.readdirSync(dir);
    for (const file of files) {
        const filePath = path.join(dir, file);
        const stat = fs.statSync(filePath);
        if (stat.isDirectory()) {
            if (!IGNORE_DIRS.includes(file)) {
                walk(filePath, fileList);
            }
        } else {
            if (!IGNORE_FILES.includes(file) && !file.endsWith('.log')) {
                fileList.push(filePath);
            }
        }
    }
    return fileList;
}

const files = walk(process.cwd());
let changedFiles = 0;
let totalReplacements = 0;

files.forEach(file => {
    try {
        let content = fs.readFileSync(file, 'utf8');
        let originalContent = content;

        SEARCH_URLS.forEach(searchUrl => {
            // Replace globally
            const regex = new RegExp(searchUrl.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'), 'g');
            content = content.replace(regex, TARGET_URL);
        });

        if (content !== originalContent) {
            fs.writeFileSync(file, content, 'utf8');
            console.log(`Updated: ${path.relative(process.cwd(), file)}`);
            changedFiles++;

            // Count replacements
            // Rough count based on length diff if needed, or just log file changed.
            // For now, simpler to just log changed files.
        }
    } catch (err) {
        console.error(`Error processing ${file}: ${err.message}`);
    }
});

console.log(`\nReplacement complete.`);
console.log(`Files modified: ${changedFiles}`);
